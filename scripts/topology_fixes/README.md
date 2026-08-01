# Topology Fixes

Scripts that convert CHARMM-GUI's PSF/PDB output into a correct GROMACS topology, applying manual repairs where CHARMM-GUI or automated tooling produced incorrect results (see `docs/KNOWN_ISSUES_AND_FIXES.md` for full diagnostic details on *why* each fix is needed).

**Which script to run depends on the system** — use the decision tree below.

## Decision tree

```
Is this a P. sativum system (WT / V53T / V83T / V53T+V83T)?
│
├── YES
│   ├── Deoxy (no O2)  →  convert_to_gmx.py            → topol_fixed.top
│   └── Oxy (with O2)  →  add_o2_leg.py                → topol_o2_leg.top
│
└── NO — O. spinosa system (WT / L43W / F125W / L43W+F125W)
    │   (these systems have the false chain-break issue — do NOT use
    │    convert_to_gmx.py's output for production; see KNOWN_ISSUES §1)
    │
    ├── Deoxy (no O2)  →  fix_chain_break.py            → topol_junction_fixed.top
    └── Oxy (with O2)  →  add_o2_and_fix_chain.py        → topol_full_fixed.top
```

## Script reference

| Script | Species / condition | Reads | Writes | Notes |
|---|---|---|---|---|
| `convert_to_gmx.py` | Any — but only used for final production on *P. sativum* systems | `step3_input.psf/.pdb` | `topol_fixed.top`, `step3_input_fixed.gro` | Plain CHARMM→GROMACS conversion via ParmEd. Do not use this output directly for *O. spinosa* systems. |
| `fix_chain_break.py` | *O. spinosa*, deoxy | `step3_input.psf/.pdb` | `topol_junction_fixed.top`, `step3_input_junction_fixed.gro` | Full conversion **plus** manual repair of the C155–N156 peptide bond/angles/dihedral. Self-contained — does not depend on `convert_to_gmx.py` having been run first. |
| `add_o2_leg.py` | *P. sativum*, oxy | `step3_input.psf/.pdb` (the base, heme-only CHARMM-GUI output — **not** a system that went through CHARMM-GUI's own O2 ligand reader) | `topol_o2_leg.top`, `step3_input_o2_leg.gro` | Builds the O2 ligand geometry and bonded parameters directly with ParmEd. |
| `add_o2_and_fix_chain.py` | *O. spinosa*, oxy | `step3_input.psf/.pdb` | `topol_full_fixed.top`, `step3_input_full_fixed.gro` | Combines the chain-junction repair and O2 construction in a single ParmEd session (avoids rebuilding the base topology twice). |

## Before running any of these

Each script expects to be run from a system's `gromacs/` working directory, with the full CHARMM36 parameter set available one level up at `../toppar/*.rtf`, `*.prm`, `*.str` (this is CHARMM-GUI's standard output layout).

```bash
cd <SYSTEM>/Solution_Builder/gromacs
python3 ../../../scripts/topology_fixes/<appropriate_script>.py
```

## Species-specific constants to check before reuse

If adapting these scripts for a new mutant/system, verify:
- `PROXIMAL_HIS_RESID` in the O2-construction scripts (93 for *P. sativum*, 110 for *O. spinosa*)
- Residue numbers `155`/`156` in the chain-break scripts (specific to the *O. spinosa* sequence used here — confirm against your own `step3_input.pdb` if reusing this for a different protein)
