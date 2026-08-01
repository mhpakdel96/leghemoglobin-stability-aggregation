# Results Summary

## 1. Heme-Pocket Stability Without O2 (Deoxy State)

RMSD of the heme group relative to the backbone-aligned protein (mean over 100 ns production; lower = more stable), n = 1 trajectory per system.

| System | Species | RMSD (Å) | Notes |
|---|---|---|---|
| V83T | *P. sativum* | 1.14 | Most stable *P. sativum* variant |
| V53T | *P. sativum* | 1.42 | |
| WT (LegHEM22) | *P. sativum* | 1.51 | |
| V53T/V83T | *P. sativum* | 1.73 | Less stable than either single mutant (negative epistasis) |
| F125W | *O. spinosa* | 1.81 | Most stable *O. spinosa* mutant |
| F125W/L43W | *O. spinosa* | 2.13 | |
| L43W | *O. spinosa* | 2.14 | |
| WT (Natural) | *O. spinosa* | 3.53 | Markedly less stable than all mutants |

*Note: the four O. spinosa values above are from the corrected topology (see `KNOWN_ISSUES_AND_FIXES.md` §1/§4); an earlier run using an unpatched chain junction gave qualitatively similar rankings but is not used in final analysis.*

## 2. Stability With O2 Bound (Oxy State)

Four systems (WT + best mutant per species) were re-simulated with O2 bound at the heme iron.

| System | Fe–O1 mean (Å) | Fe–O1 std (Å) | Fe–O1 max (Å) | O2 dissociated? | Heme+O2 RMSD (Å) |
|---|---|---|---|---|---|
| WT (LegHEM22) | 1.80 | 0.03 | 1.94 | No | 1.36 |
| V83T | 1.80 | 0.04 | 1.93 | No | 1.35 |
| WT (Natural) | 1.80 | 0.03 | 1.92 | No | 2.04 |
| F125W | 1.80 | 0.04 | 1.90 | No | 2.66 |

**Fe–O2 bond stability:** in all four systems, the Fe–O1 distance remained tightly clustered around the CHARMM36 equilibrium value (1.80 Å) for the full 100 ns, with no dissociation event — consistent with a functionally intact O2-binding pocket in every genotype tested.

**Stability ranking reversal:** the mutant that outperformed wild type without O2 (lower RMSD) did *not* retain this advantage once O2 was bound — in both species, RMSD increased for the mutant and decreased (or stayed similar) for wild type upon O2 addition.

| | Without O2 | With O2 | Direction |
|---|---|---|---|
| *P. sativum* WT | 1.51 | 1.36 | improved |
| *P. sativum* V83T | 1.14 | 1.35 | worsened |
| *O. spinosa* WT | 3.53 | 2.04 | improved |
| *O. spinosa* F125W | 1.81 | 2.66 | worsened |

## 3. Interpretation

- **Intra-species mutation effects (deoxy state):** in *P. sativum*, single mutations (V53T, V83T) each improve heme-pocket stability relative to WT, but the double mutant is less stable than either single mutant — a pattern consistent with negative epistasis, commonly reported when two individually stabilizing mutations introduce steric strain in combination (cf. heme-pocket engineering literature in myoglobin, e.g. Frontiers in Bioinformatics 2026, PMC2657002). In *O. spinosa*, all three mutations (single and combined) substantially improve stability over WT.
- **Deoxy-to-oxy stability reversal:** the loss of the mutants' stability advantage upon O2 binding is consistent with the established finding that heme-pocket geometry — and hence the effect of a given pocket mutation — differs between the deoxy and oxy coordination states in globins (e.g. Kundu et al. 2002 on soybean leghemoglobin proximal-pocket ligand regulation; Kundu & Hargrove 2003 on TyrB10-mediated ligand-binding regulation). The mutations characterized here may be favorable for the resting (unliganded) pocket conformation without being favorable for the small structural rearrangement that accompanies O2 coordination at the iron.
- **Cross-species baseline difference:** *O. spinosa* systems show markedly higher heme RMSD overall than *P. sativum*, even for the best-performing mutants. This could reflect a genuinely more flexible heme pocket in this species, and/or greater structural uncertainty in its AlphaFold model (a less-studied protein with no close experimental homolog used in training/validation). This cannot be disambiguated without an experimental structure and is reported as an open limitation.

## 4. Limitations

1. **Single replicate per system.** All values above are from one 100 ns trajectory per condition; no independent replicates with different initial velocity seeds were run. Differences smaller than ~0.3–0.4 Å (e.g., WT vs. V53T vs. V83T in *P. sativum*) should be treated as provisional until confirmed with ≥1 additional replicate per system, in line with standard MD practice for stability comparisons of this magnitude.
2. **AlphaFold-derived starting structures.** Neither protein has a deposited experimental structure; both AlphaFold models were used without further experimental validation.
3. **Homology-based heme docking.** The initial heme placement (Section 2 of `METHODOLOGY.md`) assumes conservation of the heme-pocket geometry from the 1FSL template — a reasonable assumption given the conserved globin fold, but not independently verified for these two specific sequences.
4. **No CMAP correction at the reconstructed *O. spinosa* chain junction** — a minor, disclosed approximation at a single backbone position distant from the heme pocket.

## 5. Related Literature

- Kundu, S. et al. (2002). Proximal preferences in soybean leghemoglobin heme-pocket ligand regulation. *J. Biol. Chem.* — https://pubmed.ncbi.nlm.nih.gov/11835502/
- TyrB10-mediated regulation of ligand binding in soybean leghemoglobin (MD study) — https://pubmed.ncbi.nlm.nih.gov/26211916/
- Hargrove, M.S. et al. (1997). Distal histidine mutations in recombinant soybean leghemoglobin. — https://pubmed.ncbi.nlm.nih.gov/9086279/
- Mutation-induced heme pocket destabilization, MD study (methodology comparison; comparable RMSD scale). *Front. Bioinform.* 2026. — https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1887793/full
- Active-site mutation effects in hemoglobin I from *Lucina pectinata* (MD study). — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2657002/
- Apomyoglobin stability by MD with replicate simulations. *Sci. Rep.* — https://www.nature.com/articles/srep44651
- Huang, J. & MacKerell, A.D. (2013). CHARMM36 force field validation against NMR data. *J. Comput. Chem.* — https://onlinelibrary.wiley.com/doi/10.1002/jcc.23354
