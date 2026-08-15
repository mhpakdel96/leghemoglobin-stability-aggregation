# FoldX ΔΔG Protocol

Folding free-energy calculations for all six variants.

**The FoldX binary is not redistributed here.** It is licensed software, free
for academic and non-profit use, obtained by registration at
<https://foldxsuite.crg.eu/academic-license-info>. The download link is emailed
after registration; check spam filters, as messages from the crg.es domain are
sometimes rejected by institutional filters.

Version used in this work: **FoldX 5.1, Linux 64-bit** (`foldx_20261231`,
licence valid to 2026-12-31). FoldX 5 and later no longer require
`rotabase.txt`.

---

## Sign convention

**ΔΔG = ΔG(mutant) − ΔG(wild type).**

| Sign | Meaning |
|---|---|
| Positive | Destabilising |
| Negative | Stabilising |

This is the opposite of several web-based predictors, including DynaMut2 and
the BioSig family, where positive indicates stabilisation. Mixing conventions
inverts every conclusion.

**Accuracy.** FoldX reports an intrinsic accuracy of approximately
0.5 kcal/mol. Any |ΔΔG| below that threshold is indistinguishable from zero and
should be reported as neutral, not as a small effect.

---

## Installation

```bash
mkdir -p ~/foldx5 && cd ~/foldx5
unzip ~/Downloads/foldx5_Linux_0.zip

chmod +x foldx_20261231
ln -sf "$(readlink -f foldx_20261231)" ~/foldx5/foldx
echo 'export PATH="$HOME/foldx5:$PATH"' >> ~/.bashrc
export PATH="$HOME/foldx5:$PATH"
```

Keep the `molecules/` directory alongside the binary; some commands require it.

---

## Step 1 — Verify inputs before running anything

```bash
python3 ../analysis/qc/audit_inputs.py
```

Confirm before proceeding:

| Check | Expected |
|---|---|
| Chain identifier | `A` in both models |
| *P. sativum* numbering | 1–146, residue 53 = VAL, residue 83 = VAL |
| *O. spinosa* numbering | 1–166, residue 43 = LEU, residue 125 = PHE |

A mismatch here means the mutation strings in Step 3 would silently target the
wrong residues.

---

## Step 2 — RepairPDB

Resolves minor clashes and unfavourable side-chain rotamers. **ΔΔG values
computed without this step are not valid.**

```bash
mkdir -p foldx_ddg && cd foldx_ddg
cp path/to/AF-Q9SAZ0-F1-model_v6.pdb      pea_WT.pdb
cp path/to/AF-A0A411AFI2-F1-model_v6.pdb  spi_WT.pdb

for P in pea_WT spi_WT; do
    foldx --command=RepairPDB --pdb=${P}.pdb --output-dir=.
done
```

Runtime observed: about 4 minutes per structure.

**Quality check.** RepairPDB must move side chains only, never the backbone:

| Structure | Heavy-atom RMSD | Backbone RMSD |
|---|---|---|
| pea_WT | 0.906 Å | **0.000 Å** |
| spi_WT | 1.018 Å | **0.000 Å** |

A non-zero backbone RMSD indicates a problem with the input and should be
investigated before continuing.

---

## Step 3 — Mutation files

Format: `<wild-type letter><chain><position><mutant letter>;`
Double mutants are comma-separated on a single line.

`individual_list_pea.txt`
```
VA53T;
VA83T;
VA53T,VA83T;
```

`individual_list_spi.txt`
```
LA43W;
FA125W;
LA43W,FA125W;
```

The second character is the chain identifier confirmed in Step 1.

---

## Step 4 — BuildModel

```bash
foldx --command=BuildModel --pdb=pea_WT_Repair.pdb \
      --mutant-file=individual_list_pea.txt \
      --numberOfRuns=5 --output-dir=.

foldx --command=BuildModel --pdb=spi_WT_Repair.pdb \
      --mutant-file=individual_list_spi.txt \
      --numberOfRuns=5 --output-dir=.
```

Runtime observed: 7 seconds for *P. sativum*, 7 minutes for *O. spinosa*.

**Why `--numberOfRuns=5`.** FoldX side-chain placement contains a stochastic
element. Five runs give a mean and a standard deviation instead of a single
unreproducible number.

**Why BuildModel avoids a systematic bias.** BuildModel writes both the mutant
structures and matching `WT_*` structures, and computes the difference between
them. Both therefore pass through an identical computational path. Comparing a
FoldX-generated mutant against a raw input structure would confound the mutation
effect with the effect of FoldX processing — a real and easily overlooked
artefact.

---

## Step 5 — Parse output

```bash
python3 parse_ddg.py | tee ddg_results.txt
```

Output files produced by BuildModel:

| File | Contents |
|---|---|
| `Dif_<pdb>.fxout` | Per-run energy differences, `<pdb>_<mutant>_<run>.pdb` |
| `Average_<pdb>.fxout` | Mean and standard deviation per mutant |
| `Raw_<pdb>.fxout` | Absolute energies |
| `PdbList_<pdb>.fxout` | Index of generated structures |

Row indices map to the order of lines in the mutation file: `_1_` is the first
line, `_2_` the second, `_3_` the third.

**In the `Average_*.fxout` file the column order is `Pdb`, `SD`, `total energy`
— the standard deviation precedes the mean.** Reading these in the wrong order
is an easy mistake.

---

## Results obtained

| Species | Variant | ΔΔG (kcal/mol) |
|---|---|---|
| *P. sativum* | V53T | +0.102 ± 0.029 |
| | V83T | +0.004 ± 0.0001 |
| | V53T+V83T | +0.100 ± 0.029 |
| *O. spinosa* | L43W | +1.418 ± 0.112 |
| | F125W | +0.374 ± 0.004 |
| | L43W+F125W | +1.906 ± 0.222 |

Full interpretation: [`../../docs/RESULTS.md`](../../docs/RESULTS.md).

---

## Limitations

- Calculations are on the **apo protein**, without heme. All four positions lie
  at least 8.6 Å from the cofactor, so this is a reasonable approximation, but
  it does not describe holo-protein stability.
- Input structures are AlphaFold predictions. The *O. spinosa* model
  (pLDDT 89.8) is a weaker input than the *P. sativum* model (pLDDT 95.9), and
  its ΔΔG values carry correspondingly less confidence.
- FoldX estimates folding free energy in isolation. It does not model
  aggregation, proteolysis, or expression-level effects in a living cell.

---

## Files in this directory

| File | Purpose |
|---|---|
| `individual_list_pea.txt` | Mutation definitions, *P. sativum* |
| `individual_list_spi.txt` | Mutation definitions, *O. spinosa* |
| `run_foldx.sh` | Steps 2 and 4 as a single script |
| `parse_ddg.py` | Reads `.fxout` output into a labelled table |
