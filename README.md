# Structural Stability and Aggregation Propensity of Engineered Leghemoglobin Variants

A multi-method computational assessment of solubility-directed mutations in
leghemoglobin from *Pisum sativum* and *Ononis spinosa*, combining replicated
molecular dynamics, FoldX folding free-energy calculations, and structure-based
aggregation prediction.

> **Scope correction.** This repository was originally titled *"Molecular
> Dynamics Analysis of Heme- and Oxygen-Binding Stability."* That framing was
> withdrawn after the heme cofactor was found to be harmonically restrained to
> the protein in the CHARMM-GUI topology, which makes ligand dissociation
> impossible and binding strength unmeasurable in these trajectories. Four
> further claims from the initial single-replicate analysis have also been
> retracted. See [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md).

---

## Project Question

Six site-directed mutations were designed by the originating laboratory with
two stated aims, in this order:

1. increase the stability of the mutant protein relative to wild type;
2. increase solubility to suppress inclusion-body formation during
   intracellular expression in *Pichia pastoris*.

| Species | Mutations |
|---|---|
| *P. sativum* | V53T, V83T, V53T+V83T |
| *O. spinosa* | L43W, F125W, L43W+F125W |

---

## Summary of Findings

**Stability was not increased.** In *P. sativum* all three mutations are
thermodynamically neutral (ΔΔG between +0.004 and +0.102 kcal/mol, well below
the 0.5 kcal/mol accuracy of the method). In *O. spinosa* two of the three are
genuinely destabilising (L43W +1.42, double +1.91 kcal/mol). Replicated MD
across nine independent metrics detected no significant structural effect for
any variant: of 54 mutant-versus-wild-type comparisons, two reached p < 0.05,
against ~2.7 expected by chance alone.

**Aggregation propensity was reduced — in *P. sativum* only.** All three pea
mutations lower the predicted aggregation score, additively, at essentially
zero stability cost. Residue V83 was the single highest-scoring
aggregation-prone position in the entire pea protein (+1.1027); the V83T
substitution moved it across the zero threshold to −0.7415 and shifted the
protein-wide maximum to a different residue. In *O. spinosa* only F125W helps;
L43W is detrimental, and in the double mutant the two effects cancel while the
full stability penalty remains.

**No species difference survives correction.** The originally reported
difference in heme dynamics between the two species (p = 1.3 × 10⁻⁵) was an
artefact of superimposing on a disordered N-terminal tail. Recomputed with the
structured core as the alignment reference, the difference falls from 0.60 Å to
0.04 Å (p ≈ 0.73).

### Recommended variants

| Species | Recommendation | Basis |
|---|---|---|
| *P. sativum* | **V53T+V83T**, or V83T alone | ΔA4D −8.78 (or −6.38); ΔΔG ≈ 0 |
| *O. spinosa* | **F125W only** | ΔA4D −2.10; ΔΔG +0.37. Avoid the double mutant: full +1.91 kcal/mol penalty for no aggregation benefit |

---

## Results at a Glance

| Species | Variant | ΔΔG (kcal/mol) | Δ A4D total | Verdict |
|---|---|---|---|---|
| *P. sativum* | V53T | +0.102 ± 0.029 | −2.399 | Beneficial, free |
| | V83T | +0.004 ± 0.0001 | −6.381 | Beneficial, free |
| | V53T+V83T | +0.100 ± 0.029 | **−8.780** | Best overall |
| *O. spinosa* | L43W | +1.418 ± 0.112 | +0.383 | Detrimental on both axes |
| | F125W | +0.374 ± 0.004 | −2.095 | Beneficial |
| | L43W+F125W | +1.906 ± 0.222 | −0.068 | Worst choice |

ΔΔG sign convention: positive = destabilising.
A4D sign convention: more negative = less aggregation-prone.

Full tables, statistics, and per-metric breakdowns: [`docs/RESULTS.md`](docs/RESULTS.md).

---

## Methods Overview

| Stage | Tool | Notes |
|---|---|---|
| Structure source | AlphaFold DB | Q9SAZ0 (146 aa, pLDDT 95.9); A0A411AFI2 (166 aa, pLDDT 89.8) |
| Mutagenesis | FoldX 5.1 | Backbone unchanged (RMSD 0.000); target side chain only |
| Heme placement | Superposition on 1FSL | `scripts/topology_fixes/transplant_heme.py` |
| System build | CHARMM-GUI Solution Builder | CHARMM36m, TIP3P, 0.15 M KCl, 298.15 K |
| Topology repair | ParmEd + manual patches | Three bugs fixed; see `docs/KNOWN_ISSUES.md` |
| Production MD | GROMACS 2023.4 | 13 systems × 3 replicates × 100 ns = 3.9 µs |
| Structural analysis | GROMACS + MDAnalysis | Restricted to the structured core |
| Folding stability | FoldX `BuildModel`, 5 runs | Apo protein |
| Aggregation | Aggrescan4D, pH 7.0 | Static mode, mutations applied server-side |

Detailed parameters: [`docs/METHODS.md`](docs/METHODS.md).
Web-server protocols and provenance: [`docs/WEB_TOOLS.md`](docs/WEB_TOOLS.md).

### Why analyses are restricted to a structured core

Two disordered or detached regions dominated every global metric:

- *P. sativum*: an unrepaired chain break at C143–N144 (18.30 Å) left the
  terminal tripeptide floating free. Restricting to residues 1–143 reduced
  between-replicate scatter in backbone RMSD from ±1.72 Å to ±0.06 Å — a
  28-fold improvement.
- *O. spinosa*: residues 1–13 have pLDDT < 70 and RMSF up to 18 Å. Restricting
  to residues 9–165 reduced heme-RMSD scatter three- to five-fold.

All reported structural metrics use these ranges.

---

## Repository Structure

```
├── docs/
│   ├── METHODS.md            Full parameters for every stage
│   ├── RESULTS.md            Complete result tables and statistics
│   ├── CORRECTIONS.md        Retracted claims, with evidence
│   ├── KNOWN_ISSUES.md       Seven documented bugs and artefacts
│   └── WEB_TOOLS.md          Web-server protocols and provenance
├── scripts/
│   ├── simulation/           Equilibration and production drivers
│   ├── topology_fixes/       ParmEd conversion and manual patches
│   ├── analysis/
│   │   ├── qc/               Input audit, chain integrity, flexibility
│   │   ├── structural/       Core-restricted RMSD, RMSF, Rg, SASA, H-bonds
│   │   ├── heme/             Heme RMSD, distal histidine time series
│   │   └── solubility/       Polar/apolar SASA decomposition, trend tests
│   ├── foldx/                ΔΔG protocol and parser
│   └── aggrescan4d/          Result archiving from the web server
├── results/
│   ├── md_core/              Replicated structural metrics
│   ├── sasa/                 Surface decomposition and site exposure
│   ├── heme/                 Alignment comparison, distal histidine
│   ├── foldx/                ΔΔG output
│   └── aggrescan4d/          Aggregation scores, per-residue and summary
└── structures/               AlphaFold models and repaired inputs
```

---

## Limitations

1. **Heme binding cannot be assessed.** The bonded heme model prevents
   dissociation. Metrics involving the cofactor describe the fluctuation of a
   restrained ligand.

2. **Oxygen-bound systems are not usable for ligand chemistry.** All histidines
   were built as HSD. This is correct for the proximal histidine but wrong for
   the distal histidine in the oxy form, which must be HSE to donate the
   canonical hydrogen bond to O₂. Measured consequence: that hydrogen bond is
   present in ~0 % of frames.

3. **The mutated positions lie outside the heme pocket** (8.6–12.5 Å from the
   cofactor). This was a deliberate design choice by the originating
   laboratory, not an oversight, but it means no heme-related effect was ever
   expected.

4. **Predictions, not measurements.** Neither FoldX nor Aggrescan4D has been
   validated against *Pichia pastoris* expression data. Only the pH parameter
   reflects the intended host.

5. **Single-run web analyses.** Each Aggrescan4D job was run once. The exact
   additivity observed in *P. sativum* argues for internal consistency, but
   replicate runs would provide error estimates.

6. **Three replicates limit statistical power.** For the largest observed
   effect (Cohen's d = 2.38), five replicates would be required to reach 80 %
   power.

7. **Predicted structures.** No experimental structure was available to
   validate either binding-site geometry.

---

## Reproducing This Work

Raw trajectories are not distributed here (39 files, ~336 MB each). All
analysis scripts read GROMACS output directly and are documented in
`docs/METHODS.md`. The web-server portions can be reproduced from the
checklist in `docs/WEB_TOOLS.md`.

Before any new simulation round, run the two quality-control scripts first:

```bash
python3 scripts/analysis/qc/audit_inputs.py
python3 scripts/analysis/qc/check_chain_integrity.py
```

Had these been run at the outset, the pea chain break and the mutation-position
issue would have been apparent before 3.9 µs of simulation.

---

## Software Versions

| Software | Version |
|---|---|
| GROMACS | 2023.4 |
| FoldX | 5.1 (`foldx_20261231`) |
| MDAnalysis | as installed on the analysis host |
| CHARMM force field | CHARMM36m |
| Aggrescan4D | web service, accessed 2026-08-14 |

FoldX binaries and CHARMM parameter files are licensed separately and are not
redistributed here.

---

## Licence

Code in this repository is released under the MIT Licence. Third-party
software and force-field parameters retain their own licences.
