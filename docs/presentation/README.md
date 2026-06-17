# pySeqAlign — overview presentation

`pyseqalign_overview.tex` is a Beamer/TikZ deck covering what pySeqAlign does:
relational sequence alignment (Nienhuys–Cheng distance), affine Needleman–Wunsch,
progressive multiple alignment, and relational sequence logos (ILP 2006). The
embedded `figures/*.png` are produced by the bundled examples
(`examples/scop_logos.py`, `examples/balloon_logo.py`).

Build:

```bash
latexmk -pdf pyseqalign_overview.tex      # or: pdflatex pyseqalign_overview (x2)
```

The compiled `pyseqalign_overview.pdf` is committed for convenience.
