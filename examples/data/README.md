# Example data — attribution

## SCOP (`scop/fold_*.xml`)

Protein secondary-structure-element sequences for several SCOP folds, encoded as
logical atoms `h(type,orient,len)` / `s(orient,len)` — the original data used in
Karwath & Kersting (ILP 2006). Fold files: `fold_CCC_FFF.xml` (CCC = SCOP class
1–4 → a/b/c/d, FFF = fold number); bundled folds include the TIM barrel
(`fold_003_001`), Rossmann (`fold_003_002`), immunoglobulin (`fold_002_001`),
SH3 (`fold_002_032`), and long α-hairpin (`fold_001_002`).

Derived from **SCOP** (freely available for academic use):

- Murzin A.G., Brenner S.E., Hubbard T., Chothia C. (1995). *SCOP: a structural
  classification of proteins database…* J. Mol. Biol. 247:536–540.
  <https://scop.berkeley.edu>
- Structures: **RCSB PDB** (<https://www.rcsb.org>); secondary structure via DSSP.

## Balloon (`balloon.tsv`)

Five POS-tagged sentences from Karwath & Kersting (ILP 2006), Example 1 / Table 1,
adapted from Barzilay & Lee (2003), *Learning to Paraphrase* (HLT-NAACL). Sentence 1
matches the paper's Example 1 verbatim; sentences 2–5 are tagged consistently from
Table 1's text.
