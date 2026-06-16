# pySeqAlign examples

Standalone reproductions of the **alignment + relational-logo** results from
Karwath & Kersting, *Relational Sequence Alignments and Logos* (ILP 2006).
These use only pySeqAlign (Nienhuys-Cheng atom distance + affine alignment +
multiple alignment + logos) — **no learning/boosting** is involved.

## SCOP fold logos (paper Figs 3–6)

```bash
python examples/scop_logos.py             # all bundled folds
python examples/scop_logos.py fold_003_001   # just the TIM barrel
```

Builds a multiple alignment of each fold's secondary-structure-element
sequences and renders the relational logo (stacked information-content glyphs
for the SSE atoms and their element-type abstraction, a per-column information
profile, and the per-column least-general generalisation). Output → `output/`.

## Balloon NLP example (paper Example 1 / Table 1)

```bash
python examples/balloon_logo.py
```

Aligns five POS-tagged sentences and renders the ground (`pred(word)`) and
abstract (part-of-speech) logos; total information content ≈ the paper's 80.7
bits.

## Using learned rewards (pyREAL)

The same multiple alignment accepts any scoring function, so a reward matrix
learned by pyREAL's boosting can drive the alignment instead of the fixed
Nienhuys-Cheng distance — see pyREAL's `examples/`.

Data attribution: see [`data/README.md`](data/README.md).
