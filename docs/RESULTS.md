# Results

Complete tables for all analyses. Retracted results are not reproduced here;
see [`CORRECTIONS.md`](CORRECTIONS.md).

**Systems.** 13 in total — 8 without O₂, 5 with O₂ — each run in 3 independent
100 ns replicates (39 trajectories, ≈3.9 µs aggregate). The first 20 ns of each
trajectory is discarded as equilibration.

**Analysis regions.** *P. sativum* residues 1–143; *O. spinosa* residues 9–165.
Rationale in [`METHODS.md`](METHODS.md).

---

## 1. Folding Free Energy — FoldX 5.1

`RepairPDB` followed by `BuildModel --numberOfRuns=5`. Values are the mean and
standard deviation of five runs.

**Sign convention: positive = destabilising. Method accuracy ≈ 0.5 kcal/mol;
any |ΔΔG| below that is indistinguishable from zero.**

| Species | Variant | ΔΔG (kcal/mol) | Individual runs | Verdict |
|---|---|---|---|---|
| *P. sativum* | V53T | +0.102 ± 0.029 | 0.146, 0.055, 0.110, 0.098, 0.103 | Neutral |
| | V83T | +0.004 ± 0.0001 | 0.004 ×5 | Neutral |
| | V53T+V83T | +0.100 ± 0.029 | 0.114, 0.130, 0.119, 0.088, 0.050 | Neutral |
| *O. spinosa* | L43W | +1.418 ± 0.112 | 1.630, 1.358, 1.372, 1.420, 1.310 | **Destabilising** |
| | F125W | +0.374 ± 0.004 | 0.371, 0.373, 0.380, 0.379, 0.370 | Borderline |
| | L43W+F125W | +1.906 ± 0.222 | 1.518, 2.090, 1.896, 1.874, 2.151 | **Destabilising** |

**Additivity.**

| Species | Sum of singles | Observed double | Epistasis |
|---|---|---|---|
| *P. sativum* | +0.107 | +0.100 | −0.006 |
| *O. spinosa* | +1.793 | +1.906 | +0.114 |

Both are additive within error.

**Relationship between burial and energetic cost.** Across all four positions
the ranking is strictly monotonic: the more buried the residue, the larger the
destabilisation.

| Position | Relative exposure in WT | ΔΔG |
|---|---|---|
| V83 (pea) | 70.9 % | +0.004 |
| V53 (pea) | 47.8 % | +0.102 |
| F125 (spinosa) | 37.1 % | +0.374 |
| L43 (spinosa) | 23.5 % | +1.418 |

Spearman correlation is perfect (n = 4); Pearson r = −0.84. With four points
this is illustrative rather than statistically established, but the mechanism
is straightforward — surface residues contribute little to folding stability,
so substituting them is energetically cheap.

**Limitation.** Calculations were performed on the apo protein without heme.
Since all four positions lie at least 8.6 Å from the cofactor this is a
reasonable approximation, but it does not describe holo-protein stability.

---

## 2. Aggregation Propensity — Aggrescan4D

pH 7.0 (yeast cytosol), static mode, 10 Å analysis radius, mutations applied
server-side so that wild type and mutant follow an identical path.

**Sign convention: positive residue score = aggregation-prone (highlighted by
the server); negative = solubility-promoting. A more negative total indicates
lower predicted aggregation propensity.**

| Species | Variant | Total | Average | Max | Min | Δ Total |
|---|---|---|---|---|---|---|
| *P. sativum* | WT | −139.1718 | −0.9532 | +1.1027 (V83) | −3.6158 (K115) | — |
| | V53T | −141.5708 | −0.9697 | +1.1027 (V83) | −3.6158 | −2.399 |
| | V83T | −145.5529 | −0.9969 | +1.0515 (V104) | −3.6158 | −6.381 |
| | V53T+V83T | −147.9519 | −1.0134 | +1.0515 (V104) | −3.6158 | **−8.780** |
| *O. spinosa* | WT | −173.0987 | −1.0428 | +1.9501 (V4) | −3.3658 (E100) | — |
| | L43W | −172.7161 | −1.0405 | +1.9501 (V4) | −3.3658 | +0.383 |
| | F125W | −175.1932 | −1.0554 | +1.9501 (V4) | −3.3658 | **−2.095** |
| | L43W+F125W | −173.1666 | −1.0432 | +1.9501 (V4) | −3.3658 | −0.068 |

**Additivity — the two species behave oppositely.**

| Species | Sum of singles | Observed double | Epistasis |
|---|---|---|---|
| *P. sativum* | −8.780 | −8.780 | **0.000** |
| *O. spinosa* | −1.712 | −0.068 | **+1.644** |

In *P. sativum* the two substitutions act independently and their effects sum
exactly. In *O. spinosa* they are strongly antagonistic: L43W consumes the
benefit of F125W while both stability penalties are retained.

**Site-level changes.**

| Position | WT | Mutant | Δ |
|---|---|---|---|
| V83 (pea) | V, **+1.1027** — highest score in the protein | T, −0.7415 | −1.844 |
| V53 (pea) | V, −0.4036 | T, −1.1823 | −0.779 |
| F125 (spinosa) | F, +0.3714 | W, +0.0101 | −0.361 |
| L43 (spinosa) | L, −1.0500 | W, −0.9141 | +0.136 |

V83T is the only substitution that moves a residue across the zero threshold,
from aggregation-prone to solubility-promoting. The protein-wide maximum
consequently relocates from V83 to V104.

**Aggregation-prone regions in the wild types.**

*P. sativum* — a contiguous patch spanning residues 83–106:
V83 +1.10, L84 +0.52, A91 +0.51, I92 +0.90, H93 +0.51, I94 +0.92, V99 +0.27,
V104 +1.05, V106 +0.22. Plus M1 (+0.81) at the terminus.

*O. spinosa* — dominated by the disordered N-terminal tail:
M1 +1.01, G2 +0.68, T3 +1.01, V4 +1.95, L5 +1.21. This region is invariant
across all four runs and cannot be addressed by point mutation. A secondary
cluster sits at 78–82 and 113–115.

**Cross-species comparison.** Average score per residue: *P. sativum* −0.9532
versus *O. spinosa* −1.0428, indicating that the pea protein is intrinsically
more aggregation-prone. The pea double mutant reaches −1.0134, closing most of
that gap. Note the qualitative difference: the dominant pea hotspot lies in the
structured region and is addressable; the dominant *O. spinosa* hotspot lies in
a disordered tail and is not.

---

## 2b. pH Dependence and Host Comparison

Aggrescan4D computes each job across pH 4.0–9.0. Reference points: yeast
cytosol approximately 7.0 (*Pichia pastoris*, the intended host); *E. coli*
cytosol approximately 7.6, nearest scan point 7.5.

**The Max score is blind to these mutations.** In *O. spinosa* it is identical
across all four variants at every pH point, because the highest-scoring residue
is V4 in the disordered N-terminal tail, which no mutation touches. In
*P. sativum* it shifts by at most 0.006. The server recommends Max for globular
proteins, but that recommendation does not hold when the maximum lies outside
the mutated region. Average score is used below.

### Delta versus wild type, Average score

| Variant | pH 7.0 (yeast) | pH 7.5 (*E. coli*) | Shift | Sign across pH 4–9 |
|---|---|---|---|---|
| Pea V53T+V83T | −0.0431 | −0.0420 | +0.0011 | constant |
| Pea V83T | −0.0293 | −0.0291 | +0.0002 | constant |
| Pea V53T | −0.0133 | −0.0121 | +0.0012 | constant |
| Spinosa F125W | −0.0117 | −0.0116 | +0.0001 | constant |
| Spinosa L43W | +0.0239 | +0.0265 | +0.0026 | constant |
| Spinosa double | +0.0258 | +0.0276 | +0.0018 | constant |

Additivity is reproduced independently in this module: in *P. sativum*,
−0.0133 + −0.0293 = −0.0426 against −0.0431 observed. The *O. spinosa*
antagonism is likewise reproduced: −0.0117 + +0.0239 = +0.0122 against +0.0258
observed.

### Host comparison

The benefit of each mutation is essentially pH-independent. The largest shift
between the two host conditions is 0.0026, against deltas up to 0.043. **The
same variant recommendation applies to both hosts.**

The proteins themselves, however, are predicted to be more aggregation-prone at
higher pH, both reaching a minimum near pH 6.5.

| Wild type | pH 7.0 | pH 7.5 | Change |
|---|---|---|---|
| *P. sativum* | −1.1545 | −1.1029 | +0.0516 |
| *O. spinosa* | −1.3252 | −1.2357 | +0.0895 |

**The host shift costs more than the mutations gain.** The pea double mutant at
pH 7.5 scores −1.1449, marginally worse than the wild type at pH 7.0 (−1.1545).
On this metric, moving the engineered variant from *Pichia* to *E. coli* would
forfeit the benefit obtained by mutation.

**Caveat.** This rests on the Average score. Max gives the opposite direction
for *O. spinosa* (1.8142 at pH 7.0 against 1.7944 at pH 7.5). The two metrics
disagree on host selection. Average is weighted more heavily here because Max
does not respond to the mutations at all, but this is a judgement rather than a
settled result.

Data: `results/aggrescan4d/a4d_ph_scan.csv`; raw tables in
`results/aggrescan4d/ph_scan/`.

---

## 2c. Sequence-Based Solubility Prediction (SoluProt)

An independent check using a method with no access to structure. SoluProt is
trained on soluble-expression data from *E. coli*, making it directly relevant
to the host comparison. Higher scores indicate higher predicted solubility.

Protein-Sol, the tool originally intended for this step, was unavailable on
both attempts (server error: `No space left on device`).

| Species | Variant | SoluProt | Δ vs WT |
|---|---|---|---|
| *P. sativum* | WT | 0.804 | — |
| | V53T | 0.820 | +0.016 |
| | V83T | 0.808 | +0.004 |
| | V53T+V83T | 0.817 | +0.013 |
| *O. spinosa* | WT | 0.780 | — |
| | L43W | 0.766 | −0.014 |
| | F125W | 0.794 | +0.014 |
| | L43W+F125W | 0.778 | −0.002 |

### Agreement with the structure-based analysis

The sign of every effect matches Aggrescan4D:

| Variant | SoluProt | Aggrescan4D | Agreement |
|---|---|---|---|
| Pea, all three | improved | improved | yes |
| Spinosa L43W | worse | worse | yes |
| Spinosa F125W | improved | improved | yes |
| Spinosa double | ≈ neutral | ≈ neutral | yes |

This is a meaningful independent confirmation. SoluProt uses sequence alone and
a training set drawn from bacterial expression, yet reproduces the same
qualitative conclusion — most notably the *O. spinosa* pattern in which L43W is
detrimental, F125W beneficial, and the double mutant cancels.

### Where the two methods disagree

**Ranking within *P. sativum* is close to reversed.** SoluProt places V53T
(+0.016) ahead of the double mutant (+0.013) and well ahead of V83T (+0.004);
Aggrescan4D ranks the double mutant first and V83T second, with V53T last.

The likely explanation is structural. V53T and V83T are the same substitution
(valine to threonine); a sequence-based predictor can distinguish them only
through local sequence context. The relevant difference is three-dimensional:
V83 is 70.9 % solvent-exposed and the highest-scoring aggregation hotspot in the
protein, while V53 is 47.8 % exposed and not a hotspot at all. SoluProt cannot
access this. For a target whose value derives from surface exposure, that is a
fundamental limitation, and the structure-based ranking is given more weight
here.

**Cross-species direction is opposite.** SoluProt predicts *P. sativum* (0.804)
to be more soluble than *O. spinosa* (0.780); Aggrescan4D's average score gives
the reverse. Comparing two different proteins is the least reliable use of
either tool, and no cross-species conclusion is drawn from these data.

### Resolution caveat

All eight values fall between 0.766 and 0.820, and all lie well above the 0.5
classification threshold — SoluProt predicts every variant to be soluble. The
differences between them are 0.004 to 0.016. SoluProt is a classifier rather
than a calibrated quantitative scale, and it is not established that a
difference of this size is resolvable. **The signs are reported; the precise
ranking is not.**

---

## 2d. Convergence Across Methods

Four independent lines of evidence, with different physics and different
training data, agree on the *O. spinosa* result:

| Method | L43W | F125W | Double |
|---|---|---|---|
| MD structural metrics | worse (d = 2.38 for double) | neutral | worst |
| FoldX ΔΔG | +1.42 destabilising | +0.37 borderline | +1.91 destabilising |
| Aggrescan4D | +0.383 worse | −2.095 better | −0.068 neutral |
| SoluProt | −0.014 worse | +0.014 better | −0.002 neutral |

For *P. sativum*, three methods agree that all three variants are beneficial or
neutral at no stability cost; they differ only on which single mutation
contributes most.


## 3. Structural Metrics — Replicated MD

Mean ± standard deviation across three replicates, structured core, all-atom
selections, without O₂.

| Species | Variant | RMSD (Å) | RMSF (Å) | Rg (nm) | SASA (nm²) | H-bonds | Helix (%) |
|---|---|---|---|---|---|---|---|
| *P. sativum* | WT | 1.68 ± 0.06 | 0.81 ± 0.11 | 1.51 | 88.23 ± 0.80 | 112.9 ± 0.8 | 80.97 |
| | V53T | 1.82 ± 0.21 | 0.84 ± 0.06 | 1.52 | 88.85 ± 0.40 | 112.9 ± 1.7 | 81.00 |
| | V83T | 1.70 ± 0.13 | 0.75 ± 0.05 | 1.52 | 89.56 ± 1.82 | 113.5 ± 0.9 | 79.89 |
| | V53T+V83T | 1.75 ± 0.12 | 0.86 ± 0.04 | 1.52 | 90.55 ± 0.47 | 113.4 ± 1.2 | 80.52 |
| *O. spinosa* | WT | 1.97 ± 0.14 | 0.82 ± 0.06 | 1.58 | 96.78 ± 0.81 | 123.5 ± 1.4 | 78.92 |
| | L43W | 2.02 ± 0.23 | 0.79 ± 0.05 | 1.59 | 97.36 ± 0.90 | 123.6 ± 1.3 | 79.48 |
| | F125W | 1.91 ± 0.15 | 0.85 ± 0.08 | 1.59 | 98.18 ± 0.37 | 122.2 ± 0.9 | 79.53 |
| | L43W+F125W | 2.43 ± 0.23 | 1.00 ± 0.11 | 1.59 | 98.97 ± 2.62 | 120.9 ± 1.9 | 79.27 |

Helix content of 79–81 % in every system confirms that no protein unfolded
during any trajectory.

**Statistical outcome.** Of 54 mutant-versus-wild-type comparisons (Welch's
t-test), two reached p < 0.05:

| Comparison | Δ | Cohen's d | p |
|---|---|---|---|
| Pea V53T+V83T, SASA | +2.32 nm² | 3.54 | 0.019 |
| Pea V83T, helix content | −1.08 % | −2.62 | 0.038 |

Chance alone predicts ≈2.7 results below 0.05 across 54 tests. The second
result corresponds to roughly one and a half residues of helix and is
biologically negligible. **MD alone did not resolve any mutation effect.**

Near-threshold comparisons, recorded for completeness:

| Comparison | Δ | d | p |
|---|---|---|---|
| Spinosa L43W+F125W, RMSD | +0.457 Å | 2.38 | 0.055 |
| Spinosa F125W, SASA | +1.41 nm² | 2.23 | 0.077 |
| Pea V53T, Rg | +0.007 nm | 1.96 | 0.079 |
| Spinosa L43W+F125W, RMSF | +0.181 Å | 2.02 | 0.093 |

**Species comparison after core correction** (system means, no
pseudoreplication):

| Metric | *P. sativum* | *O. spinosa* | p |
|---|---|---|---|
| RMSD | 1.74 Å | 2.08 Å | 0.056 |
| RMSF | 0.81 Å | 0.86 Å | 0.387 |
| SASA per residue | 0.625 nm² | 0.623 nm² | — |
| H-bonds per residue | 0.791 | 0.781 | — |

After normalising for chain length the two proteins are effectively
indistinguishable.

---

## 4. Surface Decomposition

Solvent-accessible surface split into apolar (C, S) and polar (N, O)
contributions. Linear regression of each metric against mutation count
(0, 1, 1, 2), 12 data points per species.

| Species | Polar slope | Apolar slope | Apolar fraction slope |
|---|---|---|---|
| *P. sativum* (V→T) | **+0.643 (p = 0.005)** | +0.521 (p = 0.093) | −0.131 (p = 0.420) |
| *O. spinosa* (→W) | +0.274 (p = 0.313) | **+0.868 (p = 0.013)** | **+0.262 (p = 0.034)** |

The two species move in opposite directions on the apolar fraction, matching
the chemistry of the substitutions: threonine adds a hydroxyl and increases
polar surface; tryptophan adds an indole ring and increases apolar surface.
This serves as an internal control confirming that the metric responds to the
intended property.

**Additivity of polar surface in *P. sativum*** is essentially exact: observed
40.882 nm² at one mutation against a linear prediction of 40.884 nm².

---

## 5. Mutation Site Exposure

Relative SASA of each target residue, normalised to the maximum accessible
surface of that amino acid, averaged over three replicates.

| Species | Position | WT residue | Relative exposure | After mutation |
|---|---|---|---|---|
| *P. sativum* | V53 | VAL | 47.8 % | THR, 42.4 % |
| | V83 | VAL | **70.9 %** | THR, 66.7 % |
| *O. spinosa* | L43 | LEU | 23.5 % | TRP, 16.6 % (buried) |
| | F125 | PHE | 37.1 % | TRP, 31.1 % |

V83 is both the most exposed target and the highest-scoring aggregation
hotspot — the ideal combination for a solubility mutation.

At position 43 the absolute SASA is unchanged (0.473 → 0.473 nm²) despite
tryptophan being substantially larger than leucine, meaning the additional
volume is accommodated inside the protein. This is consistent with the large
destabilisation measured for L43W.

---

## 6. Heme Dynamics

Heme RMSD with the structured core as alignment reference, without O₂:

| Species | WT | Mutants |
|---|---|---|
| *P. sativum* | 1.515 Å | V53T 1.479 / V83T 1.297 / double 1.717 |
| *O. spinosa* | 1.552 Å | L43W 1.453 / F125W 1.643 / double 1.500 |

No comparison is significant. Heavy-atom contact counts and Lennard-Jones
interaction energies likewise showed no effect (all p > 0.09; ΔLJ-SR between
+0.3 and +3.5 kJ/mol, all p > 0.4).

**These values describe the fluctuation of a restrained ligand and must not be
interpreted as binding strength.** See [`CORRECTIONS.md`](CORRECTIONS.md) C5.

**Heme-pocket residues (< 5 Å from the cofactor)**, provided as a basis for any
future pocket-directed design:

*P. sativum*: 32, 39, 42–45, 61, 64, 65, 68, 69, 89, 92, 93, 96, 98, 99, 102,
103, 106, 134, 137, 138, 141

*O. spinosa*: 42, 45, 46, 51–58, 75, 78, 79, 82, 102, 105, 106, 109, 110, 113,
115, 116, 119, 120, 123, 151, 154, 155, 159

Positions 61 and 93 in *P. sativum*, and 75 and 110 in *O. spinosa*, are the
distal and proximal histidines and must not be altered.

---

## 7. Distal Histidine

The only unrestrained component of the heme system, and the principal
determinant of oxygen affinity in globins.

**Effect of mutations: none.** All comparisons p > 0.091.

**Hydrogen bond to bound O₂: absent.** Across every With_O₂ system the
occupancy is approximately 0 %, with His–O₂ distances of 5.5–8.8 Å against a
canonical 2.7–3.0 Å. The cause is the HSD assignment (see
[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) item 6): with the proton on ND1, NE2
cannot donate the hydrogen bond.

**Observation, not a finding.** The *O. spinosa* distal histidine (His75) sits
outside the pocket in all 12 replicates (6.5–8.4 Å; under 1.6 % of frames below
5 Å), whereas the pea distal histidine (His61) is bistable and varies enormously
between replicates of the same system (WT: 0 %, 88.8 %, 35.6 %). The pea data
are clearly unconverged at 80 ns, and both species carry the HSD assignment, so
no conclusion is drawn.

---

## 8. Analyses That Did Not Yield Usable Results

Recorded so that they are not repeated.

**Hydrophobic patch clustering from trajectories.** Two attempts. The first
used a 10 Å centroid-to-centroid cutoff, which merged every exposed hydrophobic
residue into a single cluster spanning the whole surface. The second used
atom-to-atom contacts with a six-way parameter sweep; patch detection then
worked (2.4–8.7 clusters), but the sign of the trend inverted between exposure
thresholds in *P. sativum*, failing the stability criterion set in advance. Four
of 36 tests fell below p = 0.05 with inconsistent directions, against ≈1.8
expected by chance.

The script is retained at `scripts/analysis/solubility/patches_v2.py` with a
warning header. Aggrescan4D performs the equivalent analysis with an
experimentally calibrated scale rather than arbitrary thresholds and was used
instead.

---

*Last updated: 2026-08-15*
