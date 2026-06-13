"""Tests for Smith-Waterman local alignment."""

from pyseqalign import SmithWaterman
from pyseqalign.scoring import Blosum50, SimpleMatch


def test_identical_sequences():
    scoring = SimpleMatch(match_score=5.0, mismatch_score=-4.0, gap_score=-8.0)
    sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
    seq = [1, 2, 3, 4, 5]
    results = sw.align(seq, seq, k=1, min_score=1.0)
    assert len(results) >= 1
    assert results[0].score > 0


def test_no_match():
    scoring = SimpleMatch(match_score=5.0, mismatch_score=-4.0, gap_score=-8.0)
    sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
    seq1 = [1, 1, 1]
    seq2 = [2, 2, 2]
    results = sw.align(seq1, seq2, k=1, min_score=100.0)
    assert len(results) == 0


def test_k_zero_returns_empty():
    scoring = SimpleMatch()
    sw = SmithWaterman(scoring=scoring)
    results = sw.align([1, 2, 3], [1, 2, 3], k=0)
    assert len(results) == 0


def test_blosum50_alignment():
    """Reproduce the legacy demo: HEAGAWGHEE vs PAWHEAE."""
    scoring = Blosum50(gap_score=-8.0)
    sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
    seq1 = [9, 7, 1, 8, 1, 18, 8, 9, 7, 7]  # HEAGAWGHEE
    seq2 = [15, 1, 18, 9, 7, 1, 7]           # PAWHEAE
    results = sw.align(seq1, seq2, k=4, min_score=2.0)
    assert len(results) >= 1
    assert results[0].score > 0


def test_multiple_non_overlapping():
    scoring = SimpleMatch(match_score=10.0, mismatch_score=-4.0, gap_score=-8.0)
    sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
    # Two separate matching regions.
    seq1 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    seq2 = [1, 1, 1, 4, 4, 4, 3, 3, 3]
    results = sw.align(seq1, seq2, k=2, min_score=5.0)
    assert len(results) >= 1
