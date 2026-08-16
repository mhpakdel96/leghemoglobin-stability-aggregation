# Analysis Scripts

Organised by analysis stage. Run order is top to bottom.

> **Scripts in `legacy/` produced results that have been withdrawn.** They are
> retained only so the retracted numbers can be traced to their source. Do not
> use them for new analysis. See [`../../docs/CORRECTIONS.md`](../../docs/CORRECTIONS.md).

---

## `qc/` — quality control

**Run these before starting any simulation, not after.** Both issues they
detect were found only after 3.9 µs of production time had been committed.

| Script | Purpose |
|---|---|
| `audit_inputs.py` | Verifies model identity and organism, residue numbering, mutation correctness, and which histidine coordinates the iron |
| `check_chain_integrity.py` | Measures every C→N distance; detects spurious chain splits. A normal peptide bond is 1.33 Å |
| `diagnose_flexibility.py` | Per-residue RMSF; identifies disordered termini that must be excluded from analysis |
| `check_junction_bond.py` | Confirms the repaired PROA/PROB junction in *O. spinosa* |
| `diagnose_break.py` | Earlier chain-break diagnostic, superseded by `check_chain_integrity.py` |
| `check_by_index.py` | Atom-index lookup used during topology repair |

```bash
python3 qc/audit_inputs.py
python3 qc/check_chain_integrity.py
```

Expected: chain continuous everywhere, all four target residues confirmed as
VAL53, VAL83 (*P. sativum*) and LEU43, PHE125 (*O. spinosa*), proximal
histidine within about 2.2 Å of the iron.

---

## `structural/` — core-restricted structural metrics

Analyses are restricted to the structured core: residues 1–143 in
*P. sativum*, 9–165 in *O. spinosa*. Rationale and the 28-fold noise reduction
this produces are documented in [`../../docs/METHODS.md`](../../docs/METHODS.md)
section 4.

| Script | Purpose |
|---|---|
| `struct_analysis.sh` | Extracts protein-only trajectories; whole-chain metrics (superseded, kept for comparison) |
| `core_analysis.sh` | RMSD, RMSF, Rg, SASA, H-bonds, DSSP on the structured core |
| `fix_core_sasa_hbond.sh` | Recomputes SASA, Rg and H-bonds on all-atom core selections |
| `sasa_hydrophobic.sh` | Splits SASA into apolar (C, S) and polar (N, O) components |

**Run `fix_core_sasa_hbond.sh` after `core_analysis.sh`.** The first pass used
backbone-only selections, which exclude the side chains the mutations change
and leave `gmx hbond` with zero donors.

Submitted as SLURM jobs; edit the partition and paths at the top of each
script.

---

## `heme/` — cofactor metrics

| Script | Purpose |
|---|---|
| `hem_rmsd_core.py` | Heme RMSD under two alignment references, core and whole backbone, side by side |
| `distal_his_timeseries.py` | Distal histidine to iron distance, and to bound O₂, over the full trajectory |

`hem_rmsd_core.py` computes both alignments deliberately. The difference
between them is the origin of a retracted finding
([`../../docs/CORRECTIONS.md`](../../docs/CORRECTIONS.md) C4) and is worth
checking in any system with disordered termini.

**Interpretation limit.** Heme is harmonically restrained to the protein in
this topology. These metrics describe the fluctuation of a tethered ligand,
not binding strength.

---

## `solubility/` — surface and aggregation

| Script | Purpose |
|---|---|
| `collect_core.py` | Aggregates core metrics across all systems and replicates into CSV |
| `stats_core.py` | Welch's t-tests, mutant versus wild type, with effect sizes |
| `trend_analysis.py` | Regression of each metric against mutation count; tests additive dose response |
| `mutation_site_exposure.py` | Relative SASA of each target residue, normalised per amino acid |
| `patches_v2.py` | Hydrophobic patch clustering — **inconclusive, see header** |

`patches_v2.py` carries a warning header. Patch detection works, but the
direction of the trend inverts between exposure thresholds, failing the
stability criterion set in advance. Aggrescan4D performs the equivalent
analysis with an experimentally calibrated scale and was used instead.

---

## `legacy/` — superseded, results withdrawn

| Script | Produced | Superseded by |
|---|---|---|
| `analyze_hem_stability.py` | Single-replicate variant ranking; Fe–His "stability" | `heme/hem_rmsd_core.py` |
| `analyze_o2_stability.py` | Fe–O₁ "stability" | See `KNOWN_ISSUES.md` item 7 |

Both scripts aligned on the whole protein backbone and reported restrained
distances as results. Each carries a warning header.

---

## Typical Run Order

```bash
# 1. Before simulation
python3 qc/audit_inputs.py
python3 qc/check_chain_integrity.py

# 2. After production
sbatch --wrap="bash structural/struct_analysis.sh"
sbatch --wrap="bash structural/core_analysis.sh"
sbatch --wrap="bash structural/fix_core_sasa_hbond.sh"
sbatch --wrap="bash structural/sasa_hydrophobic.sh"

# 3. Aggregate and test
python3 solubility/collect_core.py
python3 solubility/stats_core.py
python3 solubility/trend_analysis.py
python3 solubility/mutation_site_exposure.py

# 4. Cofactor
sbatch --wrap="python3 heme/hem_rmsd_core.py"
sbatch --wrap="python3 heme/distal_his_timeseries.py"

# 5. Diagnostics, if anything looks wrong
python3 qc/diagnose_flexibility.py
```

Each script defines its system list and paths at the top; edit these to match
your directory layout.

---

## Related

| Location | Contents |
|---|---|
| [`../foldx/`](../foldx/) | ΔΔG protocol and output parser |
| [`../aggrescan4d/`](../aggrescan4d/) | Web-server result archiving |
| [`../topology_fixes/`](../topology_fixes/) | ParmEd conversion and chain repairs |
| [`../simulation/`](../simulation/) | Equilibration and production drivers |
