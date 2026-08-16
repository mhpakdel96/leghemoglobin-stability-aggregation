# Web-Based Tools: Protocols, Settings, and Provenance

This document records every web service used in this project, the exact
parameters submitted, and where the resulting data are archived in this
repository.

Web servers cannot be version-controlled. They change, go offline, or alter
their algorithms without notice. This file exists so the analyses can be
reproduced — or their limitations understood — years after the fact.

**Convention used throughout:** where a setting could not be recovered from
records, it is marked `[TO BE FILLED]` rather than guessed.

---

## 1. Structure Sources

### 1.1 AlphaFold Protein Structure Database

| Item | Value |
|---|---|
| Service | AlphaFold DB (EMBL-EBI) |
| *P. sativum* accession | `AF-Q9SAZ0-F1-model_v6` |
| *O. spinosa* accession | `AF-A0A411AFI2-F1-model_v6` |
| Downloaded | Prior to 2026-08-07 |
| Archived at | `structures/alphafold/` |

**Verification performed** (`scripts/analysis/qc/audit_inputs.py`):

| Model | UniProt annotation | Organism (TaxID) | Residues | Mean pLDDT | Residues pLDDT < 70 |
|---|---|---|---|---|---|
| Q9SAZ0 | Leghemoglobin Lb120-34 | *Pisum sativum* (3888) | 146 | 95.9 | 1 (residue 1) |
| A0A411AFI2 | Leghemoglobin | *Ononis spinosa* (58890) | 166 | 89.8 | 15 (1–13, 165–166) |

The low-confidence N-terminal segment of the *O. spinosa* model is the basis
for restricting all structural analyses to residues 9–165 in that species
(see `docs/METHODS.md`).

### 1.2 RCSB Protein Data Bank — heme placement templates

| PDB ID | Content | Use |
|---|---|---|
| `1FSL` | Ferric soybean leghemoglobin A with nicotinate | Heme transplant reference |
| `1GDJ` | Yellow lupin leghemoglobin (deoxy), 1.8 Å | Cross-check |

Heme was transplanted by structural superposition
(`scripts/topology_fixes/transplant_heme.py`), not by docking software.

**Verification:** heme centroid after transplant differs from the raw
Open Babel coordinates, confirming the superposition executed:

| File | Heme centroid (Å) |
|---|---|
| Raw `HEM.pdb` (Open Babel) | (0.17, 0.02, 0.02) |
| `LegHEM22_docked.pdb` | (0.58, 3.19, −10.12) |
| `LegHEMV53T_docked.pdb` | (0.58, 3.19, −10.12) |
| `SpinoHEMNatural_docked.pdb` | (1.74, 0.68, −8.89) |

Identical centroids across pea WT and mutants confirm a common reference
frame, which is required for a controlled comparison.

---

## 2. CHARMM-GUI — Solution Builder

The most consequential web step in the project. All simulation systems were
built here.

| Setting | Value |
|---|---|
| Module | Solution Builder |
| Force field | CHARMM36m |
| Water model | TIP3P |
| Ion concentration | 0.15 M KCl |
| Temperature | 298.15 K |
| Box type | Cubic. 6.8 nm edge for *P. sativum*, XX nm for *O. spinosa* |
| Box edge distance | approximately 1.19 nm from the protein surface |
| Ion placement method | Not recorded; CHARMM-GUI v3.7 default. 31 K⁺ and 25 Cl⁻ placed |
| Histidine protonation | **HSD for all histidines (server default)** |
| Solvent molecules | 8986 TIP3P |
| Counter-ions | 31 K⁺, 25 Cl⁻ |
| Heme parameterisation | Bonded model, standard heme patch |
| Date of build | 2026-07-28 |
| CHARMM-GUI version | 3.7 |

### Two consequences of these settings that shaped every downstream result

**Bonded heme model.** The heme iron is restrained to the four porphyrin
nitrogens (b₀ = 0.1958 nm) and to the proximal histidine NE2
(b₀ = 0.220 nm) by harmonic bonds. Heme dissociation is therefore
impossible in these trajectories. Any metric describing "heme binding
stability" measures the fluctuation of a tethered ligand, not binding
strength. See `docs/CORRECTIONS.md`.

**Histidine tautomer.** All histidines were assigned HSD (proton on ND1).
This is correct for the proximal histidine, whose NE2 must be free to
coordinate the iron. It is **incorrect for the distal histidine in the oxy
form**, which donates a hydrogen bond to bound O₂ through NE2–H. Measured
consequence: the canonical distal His–O₂ hydrogen bond is present in
approximately 0 % of frames across all With_O₂ systems, with His–O₂
distances of 5.5–8.8 Å against a canonical 2.7–3.0 Å.

Any future oxygen-binding work must rebuild the oxy systems with the distal
histidine set to HSE.

### Known CHARMM-GUI artefacts encountered

| Artefact | Detection | Fix |
|---|---|---|
| Faulty CHARMM→GROMACS conversion | Duplicate molecule counts, missing atoms | `convert_to_gmx.py` (ParmEd, explicit CHARMM36 parameters) |
| Spurious chain split, *O. spinosa* (PROA 1–155 / PROB 156–166) | C155–N156 = 1.33 Å in source model | `fix_chain_break.py` → `topol_junction_v2.top` |
| **Unrepaired chain split, *P. sativum*** | C143–N144 = 18.30 Å in WT trajectory | Not repaired; analyses restricted to residues 1–143 |

The pea chain break was discovered only after production runs were complete.
Its effect on global metrics was severe: backbone RMSD of the pea wild type
was 3.73 ± 1.72 Å over the full chain and 1.68 ± 0.06 Å over residues 1–143
— a 28-fold reduction in between-replicate scatter.

---

## 3. Aggrescan4D — aggregation propensity

The primary tool for the project's solubility question.

| Item | Value |
|---|---|
| Service | Aggrescan4D, University of Warsaw |
| URL | `https://biocomp.chem.uw.edu.pl/a4d/` |
| Date run | 2026-08-14 |
| Number of jobs | 8 (2 wild types + 6 mutants) |

### 3.1 Settings — identical for all eight jobs

| Field | Value | Rationale |
|---|---|---|
| Input | Uploaded PDB (see 3.2) | — |
| Chain(s) | A | Only chain in both models |
| pH calculations | Yes, **pH 7.0** | Yeast cytosolic pH; intracellular expression |
| Stability calculations | Yes | Provides an internal cross-check against FoldX |
| Analysis of globular regions | Yes | Excludes disordered termini |
| Dynamic mode | **No** | Static mode; dynamic jobs never left the queue |
| Distance of aggregation analysis | 10 Å | Server default |
| Improve solubility (evolutionary) | No | Automated design mode, not the question asked |
| Enhance solubility (charged) | No | Incompatible with `Mutate residues` |
| Mutate residues | Yes for mutants, **No** for wild types | See 3.3 |

### 3.2 Input files

| Species | Uploaded file | Origin |
|---|---|---|
| *P. sativum* | `pea_WT_Repair.pdb` | FoldX `RepairPDB` output |
| *O. spinosa* | `spi_WT_Repair.pdb` | FoldX `RepairPDB` output |

FoldX-repaired structures were used so that the aggregation and ΔΔG analyses
share an identical structural starting point.

### 3.3 Mutation entry — and why mutations were applied inside the server

Mutations were **not** submitted as pre-built mutant structures. The wild-type
structure was uploaded and mutations were applied by the server's
`Mutate residues` function, which calls FoldX internally.

This guarantees that wild type and mutant pass through an identical
computational path. Uploading externally built mutants would introduce a
"processed versus unprocessed" bias — the same failure mode avoided in the
FoldX `BuildModel` protocol.

Mutation strings, in the server's `<WT><new><position><chain>` notation:

| Job | Mutation string(s) |
|---|---|
| `pea_WT_pH7` | — (Mutate residues = No) |
| `pea_V53T_pH7` | `VT53A` |
| `pea_V83T_pH7` | `VT83A` |
| `pea_double_pH7` | `VT53A` + `VT83A` (single job) |
| `spi_WT_pH7` | — (Mutate residues = No) |
| `spi_L43W_pH7` | `LW43A` |
| `spi_F125W_pH7` | `FW125A` |
| `spi_double_pH7` | `LW43A` + `FW125A` (single job) |

Double mutants were submitted as **one job with two substitutions**, never as
two separate jobs.

### 3.4 Score convention

Positive residue score = aggregation-prone (highlighted yellow by the server).
Negative = solubility-promoting. A more negative total or average score
therefore indicates lower predicted aggregation propensity.

### 3.5 Results

Summary values are archived in `results/aggrescan4d/a4d_summary.csv`;
per-residue tables in `results/aggrescan4d/a4d_per_residue/`.

| Species | Variant | Total | Average | Max | Min | Δ Total vs WT |
|---|---|---|---|---|---|---|
| *P. sativum* | WT | −139.1718 | −0.9532 | +1.1027 (V83) | −3.6158 | — |
| | V53T | −141.5708 | −0.9697 | +1.1027 (V83) | −3.6158 | −2.399 |
| | V83T | −145.5529 | −0.9969 | +1.0515 (V104) | −3.6158 | −6.381 |
| | V53T+V83T | −147.9519 | −1.0134 | +1.0515 (V104) | −3.6158 | −8.780 |
| *O. spinosa* | WT | −173.0987 | −1.0428 | +1.9501 (V4) | −3.3658 | — |
| | L43W | −172.7161 | −1.0405 | +1.9501 (V4) | −3.3658 | +0.383 |
| | F125W | −175.1932 | −1.0554 | +1.9501 (V4) | −3.3658 | −2.095 |
| | L43W+F125W | −173.1666 | −1.0432 | +1.9501 (V4) | −3.3658 | −0.068 |

**Additivity.** In *P. sativum* the two single-mutant effects sum exactly to
the double-mutant effect: −2.399 + −6.381 = −8.780, matching the observed
−8.780 to three decimal places. In *O. spinosa* the effects are strongly
antagonistic: the additive prediction is −1.712 against an observed −0.068,
an epistasis of +1.644.

**Site-level changes.**

| Position | WT residue and score | Mutant residue and score | Δ |
|---|---|---|---|
| V83 (pea) | V, **+1.1027** — highest score in the protein | T, −0.7415 | −1.844 |
| V53 (pea) | V, −0.4036 | T, −1.1823 | −0.779 |
| F125 (spinosa) | F, +0.3714 | W, +0.0101 | −0.361 |
| L43 (spinosa) | L, −1.0500 | W, −0.9141 | +0.136 |

V83 was the single highest-scoring residue in the entire pea protein, and the
V83T substitution moved it across the zero threshold from aggregation-prone
to solubility-promoting. The global maximum consequently shifted from V83
(+1.1027) to V104 (+1.0515).

### 3.6 Limitations of this analysis

- Aggrescan4D has not been validated against *Pichia pastoris* expression
  data. Only the pH parameter was set to reflect the intended host.
- The predictor does not model chaperones, macromolecular crowding,
  expression rate, or degradation pathways.
- Each job was run **once**. FoldX mutant modelling contains a stochastic
  element. The exact additivity observed in *P. sativum* is strong evidence
  of internal consistency, but replicate runs would provide error estimates.
- No validated threshold exists for translating a Δ score into an expected
  change in inclusion-body formation. The only defensible threshold is the
  server's own zero line, which V83 crossed.
- The dominant hotspot in *O. spinosa* (V4, +1.9501) lies in the disordered
  N-terminal tail and is invariant across all four runs. It contributes a
  constant offset that cancels in the Δ comparison but cannot be addressed
  by point mutation.

---

## 4. Attempted Services That Did Not Yield Data

Recorded because a reader would otherwise wonder why the obvious tools are
absent, and because the failures are informative.

### 4.1 Aggrescan3D 2.0

| Item | Value |
|---|---|
| URL | `https://biocomp.chem.uw.edu.pl/A3D2/` |
| Attempted | 2026-08-14 |
| Outcome | Jobs remained `pending` for over 3.5 hours, including in static mode with dynamic mode disabled |
| Interpretation | Server-side queue stall, not a computation-time issue |

A local installation (`aggrescan3d` 1.0.2, Python 2.7 conda environment) was
attempted as a fallback. The package installed but no console entry point was
created. The effort was abandoned in favour of Aggrescan4D.

### 4.2 Protein-Sol

| Item | Value |
|---|---|
| URL | `https://protein-sol.manchester.ac.uk/` |
| Attempted | 2026-08-14 |
| Outcome | Server error: `mkdir(): No space left on device` |

Protein-Sol remains a relevant tool for a future *E. coli*-oriented analysis,
as it is trained directly on *E. coli* cell-free solubility data. It is not
appropriate as a primary tool for the *Pichia* question.

A methodological caution worth recording: sequence-based solubility
predictors rely on hydrophobicity scales in which tryptophan scores as
hydrophilic (Kyte–Doolittle −0.9) despite being among the most
aggregation-prone residues. Such tools would likely misjudge the
*O. spinosa* → W substitutions. This is precisely the limitation that
structure-based methods such as Aggrescan4D are designed to overcome.

### 4.3 FoldX licence

FoldX is not a web analysis service, but the binary is obtained through a
web registration.

| Item | Value |
|---|---|
| Registration | `https://foldxsuite.crg.eu/academic-license-info` |
| Version used | FoldX 5.1, Linux 64-bit (`foldx_20261231`) |
| Licence expiry | 2026-12-31 |

The binary is **not redistributed in this repository**. See
`scripts/foldx/README.md` for the exact commands and input files required to
reproduce the ΔΔG calculations.

---

## 5. Reproduction Checklist

To reproduce the web-based portion of this work:

1. Download both AlphaFold models using the accessions in section 1.1.
2. Run `scripts/analysis/qc/audit_inputs.py` and confirm the identity,
   numbering, and residue assignments reported there.
3. Build the simulation systems in CHARMM-GUI Solution Builder using the
   settings in section 2. **Set the distal histidine to HSE if oxygen-bound
   systems are required.**
4. Apply the topology repairs in `scripts/topology_fixes/`, then verify with
   `scripts/analysis/qc/check_chain_integrity.py`. The pea C143–N144 junction
   must be repaired in any new build.
5. Run FoldX `RepairPDB` on both wild-type models
   (`scripts/foldx/run_foldx.sh`).
6. Submit the eight Aggrescan4D jobs using the settings in section 3.1 and
   the mutation strings in section 3.3.
7. Save each per-residue table and run
   `scripts/aggrescan4d/parse_a4d_tables.py` to regenerate the summary CSV.

---

*Last updated: 2026-08-15*
