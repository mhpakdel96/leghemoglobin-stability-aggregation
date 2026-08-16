# Methods

Full parameters for every stage of the project. Web-server protocols and
provenance are documented separately in [`WEB_TOOLS.md`](WEB_TOOLS.md); known
artefacts and their mitigations in [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

---

## 1. Structure Preparation

### 1.1 Source models

Both structures are AlphaFold predictions retrieved from the AlphaFold Protein
Structure Database.

| Species | Accession | Residues | Mean pLDDT | Residues below pLDDT 70 |
|---|---|---|---|---|
| *P. sativum* | Q9SAZ0 (Leghemoglobin Lb120-34) | 146 | 95.9 | 1 |
| *O. spinosa* | A0A411AFI2 (Leghemoglobin) | 166 | 89.8 | 15 (1–13, 165–166) |

No experimental structure of either protein was available.

### 1.2 Input verification

Before any downstream work, `scripts/analysis/qc/audit_inputs.py` confirms:

- protein identity and organism from the model header;
- residue numbering continuity and the identity of every target residue;
- which histidine coordinates the heme iron;
- that each mutant differs from wild type at exactly the intended positions.

Verified assignments used throughout:

| Species | Proximal His | Distal His | Fe–NE2 (docked model) |
|---|---|---|---|
| *P. sativum* | His93 | His61 | 1.72 Å |
| *O. spinosa* | His110 | His75 | 1.26 Å |

The distal histidine of *O. spinosa* is His75, not His119. This follows from
alignment on the conserved P-K-L-x-x-H motif (pea 56–61 against *O. spinosa*
70–75) and is confirmed by proximity to the iron. The offset between the two
species is 14 residues at the distal position and 17 at the proximal position,
reflecting a three-residue insertion in *O. spinosa* between them.

### 1.3 Mutagenesis

Mutant structures were generated with FoldX 5.1 (`RepairPDB` followed by
`BuildModel`). Verified outcome: backbone RMSD between wild type and every
mutant is 0.000 Å; the only difference is the target side chain, with a maximum
atomic displacement of 2.6 Å at the substituted residue.

| Species | Variants |
|---|---|
| *P. sativum* | V53T, V83T, V53T+V83T |
| *O. spinosa* | L43W, F125W, L43W+F125W |

### 1.4 Heme placement

Heme was positioned by structural superposition onto PDB `1FSL` (ferric soybean
leghemoglobin A), not by docking software. `1GDJ` (yellow lupin leghemoglobin)
served as a cross-check.

Verification: the heme centroid is identical across pea wild type and all pea
mutants (0.58, 3.19, −10.12 Å), confirming a common reference frame. This is
required for a controlled comparison — a per-structure docking run would place
the cofactor slightly differently in each system and confound the mutation
effect.

**Note.** An earlier set of structures supplied for this project contained an
undocked heme merged at the origin, producing 79 atom pairs closer than 2 Å and
placing the iron 10.6 Å from the proximal histidine. Those files were replaced
before any simulation was run and did not enter the production pipeline.

---

## 2. System Construction

Built with CHARMM-GUI Solution Builder.

| Parameter | Value |
|---|---|
| Force field | CHARMM36m |
| Water model | TIP3P |
| Ion concentration | 0.15 M KCl |
| Temperature | 298.15 K |
| Histidine protonation | HSD for all histidines (server default) |
| Heme model | Bonded, standard heme patch |

Two default choices had consequences that persist through every result and are
documented in [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) items 6 and 7: the bonded
heme model prevents cofactor dissociation, and the HSD assignment prevents the
canonical distal His–O₂ hydrogen bond in the oxygen-bound systems.

### Topology repair

The default CHARMM-GUI output required three corrections before use, applied by
`scripts/topology_fixes/`:

1. conversion redone with ParmEd using explicit CHARMM36 parameters;
2. peptide bond restored across the spurious PROA/PROB split in *O. spinosa*;
3. a `get_bond_type` failure worked around during conversion.

A fourth, equivalent chain split in *P. sativum* (C143–N144) was not detected
until after production and was handled at the analysis stage instead.

---

## 3. Molecular Dynamics

| Parameter | Value |
|---|---|
| Engine | GROMACS 2023.4 |
| Systems | 13 (8 without O₂, 5 with O₂) |
| Replicates | 3 independent per system, different initial velocities |
| Production length | 100 ns per replicate |
| Aggregate sampling | 3.9 µs |
| Equilibration discarded | First 20 ns of each trajectory |
| Analysis frame interval | 100 ps |

Three replicates were used throughout. Single-replicate results proved
unreliable — variant rankings inverted completely when replicates were added
(see [`CORRECTIONS.md`](CORRECTIONS.md) C1).

Driver scripts: `scripts/simulation/`. The revised driver
`run_replicate_v2.sh` performs pre-flight checks; an earlier version printed a
false "Done" message on failure, silently losing 10 jobs.

---

## 4. Analysis Regions

**All structural metrics are computed on a restricted, structured core.** This
is the single most important methodological decision in the analysis.

| Species | Core (residue index) | Excluded | Reason |
|---|---|---|---|
| *P. sativum* | 1–143 | 144–146 | Detached by an unrepaired chain break (C143–N144 = 18.30 Å) |
| *O. spinosa* | 9–165 | 1–8, 166 | Disordered N-terminal tail, pLDDT below 70, RMSF up to 18 Å |

### Why this matters quantitatively

| Metric, pea wild type | Whole chain | Core only |
|---|---|---|
| Backbone RMSD | 3.73 ± 1.72 Å | 1.68 ± 0.06 Å |

Between-replicate scatter fell 28-fold. Most of the apparent variability was
the motion of a detached tripeptide.

For heme metrics the alignment reference matters as much as the measured
selection:

| Alignment region | Heme RMSD, pea | Heme RMSD, *O. spinosa* | p |
|---|---|---|---|
| Whole backbone | 1.61 Å | 2.21 Å | 1.3 × 10⁻⁵ |
| Structured core | 1.50 Å | 1.54 Å | ≈ 0.73 |

A species difference reported at p = 10⁻⁵ disappeared entirely once the mobile
tail was removed from the alignment reference. Full analysis in
[`CORRECTIONS.md`](CORRECTIONS.md) C4.

---

## 5. Structural Metrics

Computed with GROMACS tools on the structured core, all-atom selections unless
noted.

| Metric | Tool | Selection |
|---|---|---|
| Backbone RMSD | `gmx rms` | Core backbone |
| Per-residue RMSF | `gmx rmsf` | Core Cα, `-res -fit` |
| Radius of gyration | `gmx gyrate` | Core, all heavy atoms |
| Solvent-accessible surface | `gmx sasa` | Core, all heavy atoms |
| Intra-protein hydrogen bonds | `gmx hbond` | Core, all atoms including H |
| Secondary structure | `gmx dssp` | Core |

**Two selection errors worth recording**, both found and corrected: an initial
run computed SASA and Rg on backbone atoms only, which excludes the side chains
that the mutations actually change; and the hydrogen-bond group omitted
hydrogens, so `gmx hbond` reported zero donors. Both were recomputed with
all-atom core selections.

### Surface decomposition

SASA was additionally split by atom type to test the intended mechanism:

| Component | Selection |
|---|---|
| Apolar | Core atoms named `C*` or `S*` |
| Polar | Core atoms named `N*` or `O*` |

Validated to sum to the total heavy-atom count (pea: 716 apolar + 387 polar =
1103). Script: `scripts/analysis/structural/sasa_hydrophobic.sh`.

### Residue exposure

Relative SASA was normalised against maximum accessible surface values for each
amino acid (Tien et al. 2013). Classification used here: below 20 % buried,
20–50 % partially exposed, above 50 % exposed. Script:
`scripts/analysis/solubility/mutation_site_exposure.py`.

---

## 6. Statistics

**Mutant versus wild type.** Welch's t-test on replicate means (n = 3 per
system). Across 54 comparisons, two fell below p = 0.05, against approximately
2.7 expected by chance. Effect sizes are reported as Cohen's d alongside every
p-value, since with n = 3 a large effect can fail to reach significance.

**Trend tests.** Where a hypothesis predicts an additive dose response, metrics
were regressed against mutation count (0, 1, 1, 2) across all 12 replicate
values per species. This uses the full dataset and tests the hypothesis
directly, giving more power than pairwise comparisons.

**Species comparison.** Performed on system means (n = 4 per species), not on
individual replicates. Treating 12 replicate values as 12 independent samples
is pseudoreplication and inflates significance — three replicates of the same
system are not three independent draws from a species.

**Normalisation.** Cross-species comparison of extensive quantities requires
correction for chain length (143 versus 157 core residues). SASA and hydrogen
bonds are reported per residue; Rg is compared against N^0.34 scaling.

---

## 7. Folding Free Energy

FoldX 5.1. `RepairPDB` followed by `BuildModel --numberOfRuns=5`.

Sign convention: **positive ΔΔG = destabilising**. Method accuracy
approximately 0.5 kcal/mol; smaller values are reported as neutral.

`BuildModel` writes matching `WT_*` structures and computes the difference
against them, so wild type and mutant pass through an identical computational
path. Comparing a FoldX-generated mutant against a raw input structure would
confound the mutation effect with the effect of FoldX processing.

Calculations were performed on the apo protein. All four mutated positions lie
at least 8.6 Å from the cofactor, so this is a reasonable approximation, but it
does not describe holo-protein stability.

Protocol and commands: [`../scripts/foldx/README.md`](../scripts/foldx/README.md).

---

## 8. Aggregation Propensity

Aggrescan4D web server, accessed 2026-08-14.

| Setting | Value | Rationale |
|---|---|---|
| pH | 7.0 | Yeast cytosol; intracellular expression in *Pichia pastoris* |
| Mode | Static | Dynamic-mode jobs did not complete |
| Analysis distance | 10 Å | Server default |
| Globular-region analysis | Enabled | Excludes disordered termini |
| Mutations | Applied server-side | Wild type and mutant follow an identical path |
| Input structures | FoldX-repaired wild types | Shared starting point with the ΔΔG calculation |

Sign convention: **positive residue score = aggregation-prone**; negative =
solubility-promoting. A more negative total indicates lower predicted
aggregation propensity.

Neither this tool nor FoldX has been validated against *Pichia pastoris*
expression data. Only the pH parameter reflects the intended host. Full
protocol: [`WEB_TOOLS.md`](WEB_TOOLS.md) section 3.

---

## 9. Reproducibility Notes

Raw trajectories are not distributed (39 files, approximately 336 MB each). All
analysis scripts read GROMACS output directly.

**Run quality control before production, not after.** Two scripts would have
exposed the pea chain break and the mutation-position question before 3.9 µs of
simulation was committed:

```bash
python3 scripts/analysis/qc/audit_inputs.py
python3 scripts/analysis/qc/check_chain_integrity.py
```

---

## Software Versions

| Software | Version |
|---|---|
| GROMACS | 2023.4 |
| FoldX | 5.1 (`foldx_20261231`) |
| CHARMM force field | CHARMM36m |
| Aggrescan4D | Web service, accessed 2026-08-14 |
| MDAnalysis | As installed on the analysis host |

FoldX binaries and CHARMM parameter files are licensed separately and are not
redistributed in this repository.

---

*Last updated: 2026-08-16*
