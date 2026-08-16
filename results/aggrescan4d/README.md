# Aggrescan4D Results

Eight jobs: two wild types and six mutants, run at pH 7.0 in static mode.
Protocol and settings: [`../../docs/WEB_TOOLS.md`](../../docs/WEB_TOOLS.md) section 3.

| Path | Contents |
|---|---|
| `csv_export/` | Unmodified CSV exports from the web server. Filenames encode the mutation, e.g. `a3d_pea_WT [mutate_ VT53A]_A3D.csv` |
| `a4d_per_residue/` | Normalised tab-separated copies, one file per job |
| `a4d_per_residue.csv` | All jobs combined, long format |
| `a4d_summary.csv` | One row per job, recomputed from the per-residue data |
| `a4d_hotspots.csv` | Residues scoring above zero, i.e. aggregation-prone |

Summary values are **recomputed from the per-residue tables**, not copied from
the web page. Every recomputed total matched the server-reported value, which
confirms that no rows were lost during export.

**Score convention.** Positive = aggregation-prone; negative =
solubility-promoting. A more negative total means lower predicted aggregation
propensity.

Regenerate with:

```bash
python3 ../../scripts/aggrescan4d/convert_a4d_csv.py csv_export
```
