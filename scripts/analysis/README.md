# Analysis

Post-production trajectory analysis and diagnostic scripts. Run these **after** `step5_prod.tpr`/`step5_prod.xtc` exist for the systems you want to analyze.

## Script reference

| Script | Purpose | Input | Output |
|---|---|---|---|
| `analyze_hem_stability.py` | Main result: heme RMSD (backbone-aligned) for all deoxy systems | `<System>/gromacs/step5_prod.tpr/.xtc` for each system in the `SYSTEMS` list | Printed table, ranked by stability |
| `analyze_o2_stability.py` | Main result: Fe–O1 bond stability + heme+O2 RMSD for the four oxy systems | `With_O2/<System>/gromacs/step5_prod.tpr/.xtc` | Printed table |
| `check_junction_bond.py` | Diagnostic: verify the C155–N156 bond parameters actually written to a given topology file (used to confirm the fix in `topology_fixes/fix_chain_break.py` was applied correctly) | `<System>/gromacs/topol_junction_fixed.top` + `step3_input.pdb` | Printed bond-line check per system |
| `diagnose_break.py` | Diagnostic: after minimization, scan the full protein backbone for any C–N pairs farther apart than a normal peptide bond, and report which residues are affected | A `.gro` structure file (e.g. `step4.0_min.gro`) | List of broken/stretched backbone positions |
| `check_by_index.py` | Diagnostic: measure a specific bond distance using atom *serial number* rather than residue number (residue numbers are not reliable for this check after ParmEd re-numbering — see `docs/KNOWN_ISSUES_AND_FIXES.md`) | A `.pdb`/`.gro` pair | Printed distance |

## Typical usage order

1. After production finishes, run `analyze_hem_stability.py` (deoxy systems) and/or `analyze_o2_stability.py` (oxy systems) to get the headline stability numbers.
2. If a system's RMSD looks anomalously high, or you're verifying a topology fix, use `diagnose_break.py` on the minimized structure to check for any unexpectedly stretched backbone bonds, then `check_by_index.py` or `check_junction_bond.py` to pin down the exact atoms and parameters involved.

## Example

```bash
cd <path containing all system directories>
python3 /path/to/scripts/analysis/analyze_hem_stability.py
```

Edit the `SYSTEMS` list at the top of `analyze_hem_stability.py` / `analyze_o2_stability.py` to match your own directory layout and system names.
