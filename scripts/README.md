# Scripts — Overview and Execution Order

This directory contains all custom code used in the pipeline, organized by stage. Each subfolder has its own `README.md` with per-script details; this file gives the high-level order in which they run.

## Pipeline stages and where each script lives

| Stage | What happens | Scripts |
|---|---|---|
| 1. Structure prep | AlphaFold model → heme-docked PDB (PyMOL, interactive; not included as a standalone runnable script here) | *(see `docs/METHODOLOGY.md` §1–2)* |
| 2. System build | CHARMM-GUI Solution Builder (web interface) → `step3_input.psf/.pdb` | *(web-based; see `docs/METHODOLOGY.md` §3)* |
| 3. Topology conversion & repair | PSF → GROMACS topology, plus manual fixes for known CHARMM-GUI issues | [`topology_fixes/`](topology_fixes/) |
| 4. Simulation | Minimization → equilibration → 100 ns production, on GPU cluster via SLURM | [`simulation/`](simulation/) |
| 5. Analysis | Trajectory → stability metrics (RMSD, bond distances) | [`analysis/`](analysis/) |

## Overall order for a single system

```
step3_input.psf/.pdb  (from CHARMM-GUI)
        │
        ▼
[topology_fixes/ — pick the right script for this system, see topology_fixes/README.md]
        │
        ▼
topol_*.top + step3_input_*.gro
        │
        ▼
[simulation/ — SLURM job: minimization → equilibration → production]
        │
        ▼
step5_prod.tpr + step5_prod.xtc
        │
        ▼
[analysis/ — RMSD, bond-distance stability]
        │
        ▼
results/*.csv
```

Start with [`topology_fixes/README.md`](topology_fixes/README.md) — the choice of script there depends on species and whether the system includes O2, and gets the whole pipeline started correctly.
