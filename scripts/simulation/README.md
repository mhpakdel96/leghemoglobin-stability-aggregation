# Simulation (SLURM / GROMACS)

SLURM job scripts for running minimization, equilibration, and production MD on a GPU cluster.

## Script reference

| Script | Purpose | Assumes | Args |
|---|---|---|---|
| `run_full_pipeline.sh` | Generic template: minimization → equilibration → production, all in one job | Nothing pre-computed | `SYSTEM_DIR TOPOL MINGRO STARTGRO [PROD_NS]` |
| `run_production_only.sh` | Generic template: production only | `step5_prod.tpr` already generated (e.g. via a prior `gmx grompp` call) | `SYSTEM_DIR` |
| `run_o2_final.sh` | **Actual script used** for the four O2-bound systems reported in `results/`: builds a corrected `index.ndx` (via `gmx make_ndx`, merging the manually-added O2 atoms into the `SOLU` temperature-coupling group — see `docs/KNOWN_ISSUES_AND_FIXES.md` §5), then runs equilibration + production. Assumes minimization was already completed separately. | Minimized `.gro` already exists | `SYSTEM_DIR TOPOL MINGRO STARTGRO` |
| `run_spino_v2.sh` | **Actual script used** for the four corrected (chain-junction-repaired) *O. spinosa* deoxy systems: builds `index.ndx`, then runs minimization → equilibration → production in one job. | Nothing pre-computed | `SYSTEM_DIR` (paths to the `_v2` topology/coordinate files are set inside the script) |

`run_o2_final.sh` and `run_spino_v2.sh` have cluster-specific paths (GROMACS binary location, working directory root) replaced with `/path/to/...` placeholders for portability — update these for your own cluster before use.

## Example usage

```bash
# Generic full pipeline, e.g. for a new P. sativum mutant
sbatch --job-name=NewMutant scripts/simulation/run_full_pipeline.sh \
    NewMutantSystem topol_fixed.top step3_input_fixed.gro step3_input_fixed.gro 100

# Production only, once equilibration + grompp for production are already done
sbatch --job-name=NewMutant_prod scripts/simulation/run_production_only.sh NewMutantSystem
```

## Cluster resource notes

- All `mdrun` calls with GPU offload (`-nb gpu -pme gpu`) **must** run on a GPU-allocated compute node — either via `sbatch` (as here) or an interactive `srun --gres=gpu:1 --pty bash` session. Running directly on a login node fails with `Cannot run short-ranged nonbonded interactions on a GPU because no GPU is detected` (see `docs/KNOWN_ISSUES_AND_FIXES.md` §6).
- If your cluster enforces a per-user CPU/GPU quota, multiple simultaneous job submissions may queue with reason `QOSMaxCpuPerUserLimit` rather than running in parallel — this is expected and jobs will start automatically as resources free up.
