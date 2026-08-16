# Known Issues and Artefacts

Every technical problem encountered in this project, how it was detected, and
what was done about it.

Items 1–3 were found and fixed before production runs. Items 4–7 were found
afterwards and shaped how the results must be interpreted. Items 4, 5 and 7 are
the origin of the retractions in [`CORRECTIONS.md`](CORRECTIONS.md).

| # | Issue | Stage found | Status |
|---|---|---|---|
| 1 | Faulty CHARMM to GROMACS conversion | Before production | Fixed |
| 2 | Spurious chain split, *O. spinosa* | Before production | Fixed |
| 3 | `get_bond_type` failure in ParmEd | Before production | Fixed |
| 4 | Unrepaired chain break, *P. sativum* | After production | Mitigated in analysis |
| 5 | Alignment-region artefact | After production | Corrected, claim withdrawn |
| 6 | Histidine protonation state | After production | Documented, blocks O₂ work |
| 7 | Bonded heme model | After production | Documented, limits scope |

---

## 1. Faulty CHARMM to GROMACS conversion

**Symptom.** Molecule counts and atom records were inconsistent after the
default CHARMM-GUI conversion.

**Fix.** Conversion redone with ParmEd, passing the CHARMM36 parameter set
explicitly rather than relying on the server's own converter.

**Script.** `scripts/topology_fixes/convert_to_gmx.py`

---

## 2. Spurious chain split in *O. spinosa*

**Symptom.** CHARMM-GUI divided the 166-residue chain into two segments,
PROA (1–155) and PROB (156–166), and omitted the peptide bond between them.
The source model has C155–N156 = 1.33 Å, a normal peptide bond.

**Fix.** Bond, angle, and dihedral terms restored across the junction,
producing `topol_junction_v2.top`.

**Verification.** All *O. spinosa* systems, both with and without O₂, were
confirmed continuous after repair (C155–N156 = 1.30–1.34 Å in production
trajectories).

**Archived originals.** The unrepaired builds are retained on the analysis
cluster as `*_BUGGY_v1` directories. Box vectors and coordinates are identical
between `_BUGGY_v1` and `_v2` — the repair altered topology only — so any
number from the first analysis round can be traced to its source.

**Scripts.** `scripts/topology_fixes/fix_chain_break.py`,
`scripts/topology_fixes/fix_junction_batch.py`

---

## 3. `get_bond_type` failure in ParmEd

**Symptom.** Conversion aborted when resolving bond types for the heme patch.

**Fix.** Worked around within the conversion script.

**Script.** `scripts/topology_fixes/convert_to_gmx.py`

---

## 4. Unrepaired chain break in *P. sativum*

**This issue was not detected until after all production runs were complete.**

**Symptom.** The same CHARMM-GUI chain-splitting behaviour as item 2 occurred
in the pea systems, but was never repaired. The terminal tripeptide
(residues 144–146) was disconnected and drifted freely through the water box.

**Detection.** Measured C→N distance across every residue pair:

| System | C143–N144 distance |
|---|---|
| Pea wild type | **18.30 Å** |
| Pea V83T | **9.21 Å** |
| *O. spinosa* (after item 2 repair) | 1.30–1.34 Å |

A normal peptide bond is 1.33 Å. Per-residue RMSF confirmed the diagnosis:
the detached residues showed values near 30 Å, physically impossible for an
attached terminus.

**Impact.** Severe contamination of every global metric.

| Metric, pea wild type | Whole chain | Residues 1–143 |
|---|---|---|
| Backbone RMSD | 3.73 ± 1.72 Å | **1.68 ± 0.06 Å** |

Between-replicate scatter fell 28-fold. Most of the apparent variability
between replicates was the random motion of a free tripeptide, not protein
behaviour.

**Mitigation.** All structural analyses are restricted to residues 1–143. The
artefact is identical across every pea system, so mutant-versus-wild-type
comparisons within *P. sativum* remain internally consistent.

**Why systems were not rebuilt.** The break lies at the extreme C-terminus,
two residues beyond the last heme-contacting residue (141), and outside the
region under study. Rebuilding six systems and rerunning 18 trajectories was
judged not worth the marginal gain. This is a documented compromise, not an
oversight.

**Caveat.** Residues 144–146 have pLDDT above 70 in the AlphaFold model and
lie packed against the protein surface (23 of 32 heavy atoms within 4 Å of the
core). Detaching them is a local, peripheral perturbation — not a harmless one.

**Detection script.** `scripts/analysis/qc/check_chain_integrity.py`

---

## 5. Alignment-region artefact

**The most consequential issue in the project.** Full analysis in
[`CORRECTIONS.md`](CORRECTIONS.md) C4.

**Symptom.** Heme RMSD appeared dramatically higher in *O. spinosa* than in
*P. sativum*, with p = 1.3 × 10⁻⁵ across two metrics and two statistical tests.

**Cause.** RMSD of any molecule requires prior superposition on a reference.
The analysis superimposed on the entire protein backbone, which in
*O. spinosa* includes residues 1–8 with pLDDT below 70 and RMSF up to 18 Å.
Aligning on a mobile reference attributes spurious motion to everything else.

| Alignment region | *P. sativum* | *O. spinosa* | Difference | p |
|---|---|---|---|---|
| Whole backbone | 1.61 Å | 2.21 Å | 0.60 Å | 1.3 × 10⁻⁵ |
| Structured core | 1.50 Å | 1.54 Å | 0.04 Å | ≈ 0.73 |

**Independent confirmation.** Between-replicate scatter in *O. spinosa* fell
from ±0.35 Å to ±0.07 Å under core alignment. A five-fold noise reduction is
what is expected when a noise source is removed, not when a real signal is
suppressed.

**Mitigation.** All heme metrics use the structured core as alignment
reference: residues 1–143 in *P. sativum*, 9–165 in *O. spinosa*.

**General lesson.** In any system containing disordered termini, the choice of
alignment region should be treated as an analysis parameter and tested, not
assumed. It can manufacture a species-level difference with a p-value of 10⁻⁵.

**Scripts.** `scripts/analysis/heme/hem_rmsd_core.py` computes both alignments
side by side; `scripts/analysis/qc/diagnose_flexibility.py` identifies
disordered regions.

---

## 6. Histidine protonation state

**Symptom.** The canonical hydrogen bond between the distal histidine and
bound O₂ is present in approximately 0 % of frames across all With_O₂ systems.
His–O₂ distances range from 5.5 to 8.8 Å, against a canonical 2.7–3.0 Å.

**Cause.** CHARMM-GUI assigned HSD (proton on ND1) to every histidine. This is
correct for the proximal histidine, whose NE2 must remain free to coordinate
the iron. It is wrong for the distal histidine in the oxy form, which donates
a hydrogen bond to bound O₂ through NE2–H. With the proton on ND1, NE2 has no
hydrogen to donate, and the two groups repel.

**Impact.** The oxygen-bound systems do not reproduce the principal
protein–ligand interaction of an oxy-globin. **No conclusion about oxygen
binding or affinity can be drawn from them.**

**Fix required.** Rebuild all five With_O₂ systems in CHARMM-GUI with the
distal histidine set to HSE — residue 61 in *P. sativum*, residue 75 in
*O. spinosa*. This is a single selection in the protonation-state step of the
build, but it requires new production runs.

**Note on residue identity.** The distal histidine in *O. spinosa* is His75,
not His119. Established by sequence alignment on the conserved P-K-L-x-x-H
motif (pea 56–61 aligns with *O. spinosa* 70–75) and confirmed by proximity to
the heme iron.

**Script.** `scripts/analysis/heme/distal_his_timeseries.py`

---

## 7. Bonded heme model

**Not a bug.** This is standard, correct CHARMM-GUI practice for heme
proteins. It is recorded here because it defines the boundary of what these
simulations can measure.

**Description.** The heme iron is harmonically restrained to the four
porphyrin nitrogens (b₀ = 0.1958 nm) and to the proximal histidine NE2
(b₀ = 0.220 nm). In the oxy systems, Fe–O₁ is likewise restrained
(b₀ = 0.180 nm).

**Consequence.** Heme dissociation and oxygen release are impossible in these
trajectories, for any variant, under any analysis.

**How this was detected.** Observed distances carried a standard deviation of
**exactly zero** across all systems and replicates. No physical quantity in a
finite-temperature simulation has zero variance. The values were the input read
back out.

```
[ bonds ]
2289  2216   1   0.18000  209200.0     ; Fe-O1,        b0 = 0.180 nm
1430  2218   1   0.22000   54392.0     ; Fe-NE2(His),  b0 = 0.220 nm
```

**Impact.** Two reported "findings" were withdrawn on this basis
([`CORRECTIONS.md`](CORRECTIONS.md) C2 and C3), and the scope of the project
was revised (C5).

**Alternatives, if heme binding must be assessed.** Reparameterise heme
non-bonded and apply MM-PBSA or an alchemical free-energy method. Neither was
pursued, as all four mutated positions lie 8.6–12.5 Å from the cofactor and no
heme-related effect was expected.

---

## Legacy Scripts

Two analysis scripts in `scripts/analysis/legacy/` produced results that have
since been withdrawn. They are retained so that the retracted numbers can be
traced to their source, and carry warning headers.

| Script | Produced | Status |
|---|---|---|
| `analyze_hem_stability.py` | Single-replicate variant ranking; Fe–His "stability" | Superseded by `heme/hem_rmsd_core.py` |
| `analyze_o2_stability.py` | Fe–O₁ "stability" | Superseded; see item 7 |

Neither should be used for new analysis.

---

## Process Recommendation

Two quality-control scripts should be run **before** any future production
round, not after:

```bash
python3 scripts/analysis/qc/audit_inputs.py
python3 scripts/analysis/qc/check_chain_integrity.py
```

The first verifies model identity, residue numbering, mutation correctness, and
the proximal-histidine assignment. The second detects chain breaks.

Run at the outset, they would have exposed items 4 and 6 before 3.9 µs of
simulation time was committed.

---

*Last updated: 2026-08-16*
