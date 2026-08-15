# Corrections and Retractions

This document records claims made in earlier versions of this project that
were subsequently found to be unsupported, together with the evidence that
overturned them.

It is kept deliberately visible. Four of the five entries below arose from
technical artefacts that are common in molecular dynamics of cofactor-bound
proteins and are rarely reported. Documenting them is more useful than
quietly rewriting the history.

**Status of all entries: withdrawn. None of the claims below should be cited.**

---

## C1 — Ranking of mutant variants by heme stability

**Withdrawn claim.** In *O. spinosa*, the double mutant F125W+L43W stabilised
heme binding most effectively, and the wild type was the least stable variant.
A "mirror-image pattern" was reported between the two species, with single
mutations improving stability in *P. sativum* and the combination worsening it.

**Basis of the claim.** One 100 ns replicate per system.

**Evidence against.** With three independent replicates per system, the
ranking inverted:

| System | Replicate 1 only | Mean of 3 replicates |
|---|---|---|
| WT | 2.03 Å (rank 2) | 2.34 Å (worst) |
| L43W | 2.40 Å (rank 3) | 1.99 Å (best) |
| F125W | 2.43 Å (worst) | 2.21 Å (rank 3) |
| F125W+L43W | 1.75 Å (best) | 2.08 Å (rank 2) |

The best-ranked variant became second, the worst became best, and the wild
type moved from second to last. The differences were within replicate noise.

**Lesson.** A single trajectory is a single draw from a distribution. Random
variation almost always produces an apparent pattern.

---

## C2 — Stability of the Fe–O₁ bond

**Withdrawn claim.** The Fe–O₁ distance remained stable at 1.80 Å across all
oxygen-bound systems, indicating that the bound oxygen never dissociated.

**Evidence against.** The value is a harmonic restraint written into the
topology:

```
[ bonds ]
2289  2216   1   0.18000  209200.0     ; Fe–O1, b0 = 0.180 nm
```

The observed standard deviation was **exactly zero** across all systems and
replicates. No physical quantity in a finite-temperature MD simulation has
zero variance. The number is the input read back out.

---

## C3 — Stability of the Fe–His distance

**Withdrawn claim.** The Fe–NE2(His) distance remained at 2.24 ± 0.06 Å in
all 13 systems, demonstrating a stable proximal coordination bond.

**Evidence against.** Same mechanism as C2:

```
1430  2218   1   0.22000   54392.0     ; Fe–NE2(His), b0 = 0.220 nm
```

The observation confirms only that the restraint was applied and that
minimisation resolved the initial clash — the docked structures began at
1.72 Å (pea) and 1.26 Å (*O. spinosa*), both physically impossible.

---

## C4 — Species difference in heme dynamics

**This is the most consequential retraction.**

**Withdrawn claim.** Heme is significantly more mobile in *O. spinosa*
leghemoglobin than in *P. sativum* (2.15 Å vs 1.55 Å, p = 1.3 × 10⁻⁵),
confirmed by two independent metrics and two statistical tests.

**Evidence against.** Computing the RMSD of any molecule requires first
superimposing on a reference. The original analysis superimposed on the
entire protein backbone. The *O. spinosa* backbone includes an N-terminal
segment (residues 1–8) with pLDDT below 70 and RMSF up to 18 Å. Aligning on
a mobile reference attributes spurious motion to everything else.

| Alignment region | *P. sativum* | *O. spinosa* | Difference | p |
|---|---|---|---|---|
| Whole backbone | 1.61 Å | 2.21 Å | 0.60 Å | 1.3 × 10⁻⁵ |
| Structured core only | **1.50 Å** | **1.54 Å** | **0.04 Å** | **≈ 0.73** |

Ninety-four per cent of the reported difference disappeared. Independent
confirmation: between-replicate scatter in *O. spinosa* fell from ±0.35 Å to
±0.07 Å under core alignment — a five-fold reduction, exactly what is
expected when a noise source is removed rather than a signal.

The protein-level species difference likewise did not survive. After
restricting to the structured core and normalising for chain length
(143 versus 157 residues), per-residue SASA (0.625 vs 0.623 nm²) and hydrogen
bonds per residue (0.791 vs 0.781) are effectively identical, and the RMSF
difference is not significant (p = 0.387).

**Lesson.** The choice of alignment region can manufacture a species-level
difference with a p-value of 10⁻⁵. This failure mode is specific to systems
containing disordered termini and is worth checking routinely.

---

## C5 — Project scope: "heme- and oxygen-binding stability"

**Withdrawn framing.** The repository was originally presented as measuring
the stability of heme and oxygen binding.

**Evidence against.** CHARMM-GUI builds heme proteins with the iron
harmonically restrained to the four porphyrin nitrogens (b₀ = 0.1958 nm) and
to the proximal histidine NE2 (b₀ = 0.220 nm). This is standard, correct
practice for keeping a cofactor in place — but it makes dissociation
impossible. Binding strength therefore cannot be measured in these
trajectories, for any variant, under any analysis.

Additionally, all four mutated positions lie 8.6–12.5 Å from the cofactor,
outside the heme pocket. Van der Waals interactions decay as r⁻⁶, so no
direct effect on heme binding was physically expected. This was a deliberate
design decision by the originating laboratory — the mutations were selected to
avoid the catalytic site — and is not itself an error.

**Corrected scope.** The project assesses structural stability of the folded
state, folding free energy, and predicted aggregation propensity. See
[`README.md`](../README.md).

---

## What Was Verified and Found Correct

A full audit of the input chain was performed
(`scripts/analysis/qc/audit_inputs.py`). The following were confirmed:

| Item | Result |
|---|---|
| Protein identity | Both models are annotated leghemoglobins from the intended species |
| Residue numbering | No offset. His93 is proximal (Fe–NE2 1.72 Å in the docked model); His61 distal |
| Mutation correctness | Each mutant differs from wild type at exactly the intended position(s), with no unintended substitutions |
| Docking templates | 1FSL and 1GDJ are both genuine leghemoglobin structures |
| Heme placement | Superposition executed correctly; identical heme frame across wild type and mutants |
| Wild-type/mutant symmetry | Backbone RMSD 0.000 between wild type and every mutant; only the target side chain differs |
| CHARMM-GUI input | Rigid-body transform only (RMSD 0.000 after superposition); wild type and mutant systems built identically |
| *O. spinosa* topology | Chain continuous after repair (C155–N156 = 1.30–1.34 Å) |

An earlier concern that the colleague-supplied structures contained an
undocked heme (79 atom pairs closer than 2 Å; Fe 10.6 Å from the proximal
histidine) was confirmed, but those files were superseded by a repeated
docking procedure before any simulation was run. They did not enter the
production pipeline.

---

*Last updated: 2026-08-15*
