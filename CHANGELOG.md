# Changelog

## 0.1.0

Initial release -- pure Python rewrite of the legacy pyAlign and pySeqAlign libraries.

### Alignment
- Smith-Waterman local alignment with k-best non-overlapping results
- Needleman-Wunsch global alignment
- BLOSUM50 substitution matrix
- Simple identity-match scorer
- Nienhuys-Cheng atom distance for structured atoms
- Optional SWI-Prolog integration via pyswip
- Amino acid encoding/decoding utilities
- Bundled Prolog knowledge base files (distances, BLOSUM50, amino acids)

### Learning (ILP)
- Common ILP task/result types and learner protocol
- AlignmentTaskBuilder for converting sequence pairs to ILP format
- Aleph backend (classic ILP, ported from legacy pySeqAlign)
- Popper backend (modern ILP, learning from failures)
- Dual file output: Aleph (.b/.f/.n) and Popper (bk.pl/exs.pl/bias.pl)
- Bundled aleph_swi_ak.pl for SWI-Prolog compatibility
