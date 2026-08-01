# Known Issues Encountered During Pipeline Development, and Their Resolutions

This document records software/tooling issues discovered during pipeline development, the diagnostic evidence for each, and the fix applied. It is included for transparency and reproducibility — several of these are not documented elsewhere and may be useful to others building similar pipelines.

## 1. False chain-break detection by CHARMM-GUI (*O. spinosa* systems)

**Symptom:** CHARMM-GUI's Model/Chain Selection page split the *O. spinosa* protein into two protein segments (PROA: residues 1–155, PROB: residues 156–166) instead of one continuous chain.

**Root-cause investigation:** direct measurement of the C(155)–N(156) bond distance in (a) the raw, unprocessed AlphaFold model, (b) the pre-CHARMM-GUI heme-docked PDB, both gave 1.33 Å — identical to backbone bond lengths measured elsewhere in the same chain (e.g., C50–N51, C100–N101, C160–N161, all 1.34 Å) and to the accepted peptide bond length. No `TER` record or chain-ID change was present at this position in the docked PDB fed to CHARMM-GUI. **Conclusion: the input coordinates and connectivity were correct; the split is an artifact of CHARMM-GUI's own segment-assignment heuristic**, not a defect in the AlphaFold model or the docking script.

**Fix:** Terminal patching for both PROA and PROB was set to `NONE` during CHARMM-GUI processing (preventing incorrect N-/C-terminal charge patches at the artificial junction), and the missing peptide bond, associated angles, and omega dihedral were reconstructed post hoc with ParmEd using explicit CHARMM36 parameter values (see `METHODOLOGY.md`, §4b).

## 2. CHARMM-GUI auto-generated force field for O2 ligand read as a separate hetero chain

**Symptom:** when the O2 ligand was submitted to CHARMM-GUI as a structurally distinct "Hetero" chain (rather than embedded within the same segment as the heme group), CHARMM-GUI's ligand reader did not match it against the existing `RESI O2` definition in `toppar_all36_prot_heme.str`, and instead auto-generated a new force field for it. The resulting single-point energy was anomalously high (~7.5×10⁷ kJ/mol VDW term, vs. <2×10⁷ for all other systems), indicating a severe structural/parametrization problem.

**Fix:** abandoned the CHARMM-GUI ligand-reader pathway for O2 entirely. O2 atoms are instead constructed directly in the GROMACS-format topology using ParmEd, after the base system (protein + heme, no O2) has already been correctly built via CHARMM-GUI (see `METHODOLOGY.md`, §4c).

## 3. Duplicate protein molecule count in CHARMM-GUI's native GROMACS converter

**Symptom:** CHARMM-GUI's built-in CHARMM→GROMACS conversion produced a `topol.top` listing the protein molecule with multiplicity 2 (`PROA 2`) despite only one copy existing, and with an atom count mismatch versus the PSF (2290 vs. 2322 atoms).

**Fix:** bypassed CHARMM-GUI's converter; topology generated directly from the PSF with ParmEd (`convert_to_gmx.py`).

## 4. Silent incorrect parameter lookup (`get_bond_type` / `get_angle_type` helper functions)

**Symptom:** a custom Python helper function that queried `CharmmParameterSet.bond_types[(type1, type2)]` intermittently returned an internally consistent-looking but **numerically incorrect** `BondType`/`AngleType` object — not a `KeyError`, and not the requested parameter. This was discovered independently for two unrelated bonds (Fe–O1 in the O2 ligand, and the C(155)–N(156) peptide bond above) within the same Python session; in both cases the function returned identical, wrong values (b0 = 0.180 nm, k = 209200 kJ·mol⁻¹·nm⁻², which correspond to neither bond's correct parameters).

**Diagnostic method:** cross-checking the *n*-th bond's parameters as written in the final `.top` file against the raw CHARMM parameter file value for the correct atom-type pair, using the atom serial numbers (not residue numbers, which can be renumbered during PSF→GROMACS conversion).

**Fix:** all custom-added bonded terms (bonds, angles, dihedrals) in this pipeline are now parametrized with explicit, hand-verified numeric values read directly from the raw CHARMM36 `.prm`/`.str` files, rather than via automated `ParameterSet` dictionary lookups. This is the standard used throughout `add_o2_and_fix_chain.py` and `add_o2_leg.py`.

## 5. GROMACS temperature-coupling group omission after manual atom insertion

**Symptom:** `gmx grompp` failed with `Fatal error: 2 atoms are not part of any of the T-Coupling groups` when equilibrating systems where O2 atoms had been added directly via ParmEd (bypassing CHARMM-GUI's own index-file generation, which therefore did not know about the new atoms).

**Fix:** regenerate `index.ndx` per system with `gmx make_ndx`, explicitly merging the new atoms into the appropriate temperature-coupling group (`Protein | HEME | O2` → `SOLU`; `POT | CLA | TIP3` → `SOLV`).

## 6. `gmx mdrun -nb gpu` fails with "no GPU is detected" when run on a SLURM login node

**Symptom:** `Fatal error: Cannot run short-ranged nonbonded interactions on a GPU because no GPU is detected` when `mdrun` was invoked interactively without first requesting a GPU-allocated compute node.

**Fix:** GPU-accelerated `mdrun` steps must be run either inside an `srun --gres=gpu:a100:1 --pty bash` interactive session on a compute node, or submitted as an `sbatch` job requesting a GPU resource — never directly on the cluster login node.
