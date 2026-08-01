# Methodology

Complete, reproducible protocol for all 16 simulated systems (8 variants × ±O2).

## 1. Starting Structures

- *P. sativum* leghemoglobin: AlphaFold2 model `AF-Q9SAZ0-F1-model_v6.pdb` (UniProt Q9SAZ0).
- *O. spinosa* leghemoglobin: AlphaFold2 model `AF-A0A411AFI2-F1-model_v6.pdb` (UniProt A0A411AFI2).

No experimental (crystallographic) structure was available for either protein; this is disclosed as a limitation (see `RESULTS_SUMMARY.md`).

Point mutants were generated in silico (FoldX repair/mutate protocol) prior to heme docking.

## 2. Heme Cofactor Docking

**Tool:** PyMOL (script: `scripts/topology_fixes/transplant_heme.py` — not automated in this repo, run interactively).

**Method:** homology-based cofactor transplantation. PDB entry `1FSL` (a leghemoglobin structure with heme bound) was used as template. For each target chain, the pipeline:
1. Superimposes the target chain onto template chain A/B (whichever gives shorter Fe–proximal-His distance).
2. Extracts the HEM heteroatom group from the template.
3. Merges it into the target structure at the corresponding pocket position.

**Output:** `<System>_docked.pdb` — protein + heme, Fe positioned near (but not yet covalently patched to) the proximal histidine. Initial Fe–His distances (1.3–1.7 Å before minimization) reflect only geometric transplantation, not the true equilibrium bond length; this is corrected during minimization (Section 5) and by the PHEM patch (Section 3).

## 3. System Preparation (CHARMM-GUI Solution Builder)

Web interface: [charmm-gui.org](https://charmm-gui.org) → Input Generator → Solution Builder.

**Steps performed for every system:**
1. **PDB Reader & Manipulator** — upload `<System>_docked.pdb`; apply the `PHEM` patch (heme–proximal-His coordination) with the correct His residue number (93 for *P. sativum*, 110 for *O. spinosa* — CHARMM-GUI's automatic guess was always incorrect and required manual override).
2. **Solvator** — rectangular box, 10 Å edge distance, TIP3P water, 0.15 M KCl (Monte Carlo ion placement).
3. **PBC Setup** — periodic boundary conditions, automatic PME grid.
4. **Input Generator** — CHARMM36m force field, GROMACS output format, NVT equilibration + NPT production ensembles, target temperature 298.15 K.

**Key outputs used downstream:** `step3_input.psf`, `step3_input.pdb`, `step4.0_minimization.mdp`, `step4.1_equilibration.mdp`, `step5_production.mdp`, and the full CHARMM36 parameter set (`toppar/`).

### 3a. *O. spinosa*-specific chain handling

CHARMM-GUI's chain-detection heuristic erroneously split the *O. spinosa* protein into two segments (PROA: residues 1–155, PROB: residues 156–166) despite the source PDB showing a fully continuous, correctly bonded backbone at this junction (verified directly against the raw AlphaFold coordinates and the pre-CHARMM-GUI docked file — see `docs/KNOWN_ISSUES_AND_FIXES.md`). Terminal patching for both segments was set to `NONE` at this stage to avoid CHARMM-GUI erroneously applying N-/C-terminal charge patches at the artificial junction; the correct peptide bond was reconstructed later (Section 4b).

## 4. Topology Conversion and Manual Repair

**Tool:** ParmEd (Python), using the full CHARMM36 parameter set from `toppar/`.

### 4a. Base conversion (`scripts/topology_fixes/convert_to_gmx.py`)

CHARMM-GUI's built-in CHARMM→GROMACS converter produced an incorrect topology (duplicate protein molecule count, missing atoms). ParmEd was used instead to load the PSF directly with the full parameter set and export a correct GROMACS `.top`/`.gro` pair.

### 4b. Chain-junction repair (*O. spinosa* systems only)

**Script:** `scripts/topology_fixes/fix_chain_break.py` (initial version) / superseded by the manually-parametrized version integrated into `add_o2_and_fix_chain.py`.

Adds the missing C(155)–N(156) peptide bond plus associated angles (CA-C-N, O-C-N, C-N-CA, C-N-H) and the omega backbone dihedral, using values read directly from `par_all36m_prot.prm`:

| Term | Parameters (CHARMM units) |
|---|---|
| Bond C–N | Kb = 370.0 kcal·mol⁻¹·Å⁻², b0 = 1.345 Å |
| Angle CA-C-N | K = 80.0, θ0 = 116.5° |
| Angle O-C-N | K = 80.0, θ0 = 122.5° |
| Angle C-N-CA | K = 50.0, θ0 = 120.0° |
| Angle C-N-H | K = 34.0, θ0 = 123.0° |
| Dihedral ω (CA-C-N-CA) | Two terms: (K=1.6, n=1, φ=0°), (K=2.5, n=2, φ=180°) |

**Important:** an automated CHARMM-parameter-set lookup function (`get_bond_type`/`get_angle_type`) was found to intermittently return incorrect parameter objects (see `KNOWN_ISSUES_AND_FIXES.md`). All bonded terms added by custom scripts in this pipeline therefore use explicit, manually verified numeric values rather than automated lookups.

### 4c. O2 ligand construction (`scripts/topology_fixes/add_o2_leg.py`, `add_o2_and_fix_chain.py`)

The O2 ligand is **not** read through CHARMM-GUI's ligand reader (which mis-parametrized it via automatic force-field generation). Instead, O2 is constructed directly with ParmEd on the existing (heme-only) topology:

1. Position O1 along the vector from Fe opposite the proximal His, at 1.80 Å (equilibrium Fe–O bond length, `toppar_all36_prot_heme.str`).
2. Position O2 at a 120° Fe–O1–O2 bend, O1–O2 = 1.23 Å.
3. Add atoms with CHARMM atom type `OM` (heme-ligand oxygen), charges +0.021/−0.021.
4. Add bonded terms with explicit parameters:

| Term | Parameters |
|---|---|
| Bond O1–O2 | Kb = 600.0, b0 = 1.23 Å |
| Bond Fe–O1 | Kb = 250.0, b0 = 1.80 Å |
| Angle O2-O1-Fe | K = 0 (intentionally zero per CHARMM36 heme parameter file), θ0 = 180° |
| Angle O1-Fe-N(porphyrin) ×4 | K = 5.0, θ0 = 90° |
| Dihedral O2-O1-Fe-N | K = 0 (intentionally zero) |

For *O. spinosa* systems, the chain-junction repair (4b) and O2 construction (4c) are performed together in a single ParmEd session (`add_o2_and_fix_chain.py`) to avoid regenerating the base topology twice.

## 5. Minimization

GROMACS `mdrun`, steepest-descent, Fmax < 1000 kJ·mol⁻¹·nm⁻¹ convergence criterion, ≤5000 steps.

## 6. Equilibration

NVT ensemble, 125 ps, 298.15 K, Maxwell-Boltzmann-distributed initial velocities (random seed, `gen_seed = -1` in GROMACS `.mdp`).

**Note on index groups:** for systems where O2 was added by direct ParmEd atom insertion, the original CHARMM-GUI-generated `index.ndx` did not include the new atoms in any temperature-coupling group. A new index file was generated per system with `gmx make_ndx`, merging Protein+HEME+O2 into a `SOLU` group and POT+CLA+TIP3 into `SOLV`.

## 7. Production

NPT ensemble, 100 ns (50,000,000 steps × 2 fs), 298.15 K, GROMACS 2023.4 on NVIDIA A100 GPUs (SLURM cluster). Representative job scripts: `scripts/simulation/`.

## 8. Trajectory Analysis

**Tool:** MDAnalysis (Python).

- **Heme/O2–protein RMSD:** for each frame, the protein backbone is aligned to the first production frame (translation + rotation only); RMSD of the heme (or heme+O2) atom positions is computed against their first-frame positions in this aligned reference frame. This measures cofactor displacement relative to the protein scaffold, independent of overall protein tumbling in the simulation box.
- **Fe–His / Fe–O1 bond distance:** simple Euclidean distance between the specified atom pair, computed per frame across the full trajectory.

Analysis scripts: `scripts/analysis/`.
