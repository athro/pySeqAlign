"""Tests for Needleman-Wunsch global alignment."""

from pyseqalign import NeedlemanWunsch
from pyseqalign.scoring import Blosum50, SimpleMatch


def test_identical_sequences():
    scoring = SimpleMatch(match_score=5.0, mismatch_score=-4.0, gap_score=-8.0)
    nw = NeedlemanWunsch(scoring=scoring)
    seq = [1, 2, 3, 4, 5]
    result = nw.align(seq, seq)
    assert result.query == seq
    assert result.target == seq
    assert result.score > 0


def test_gap_insertion():
    scoring = SimpleMatch(match_score=5.0, mismatch_score=-4.0, gap_score=-2.0)
    nw = NeedlemanWunsch(scoring=scoring, gap_penalty=-2.0)
    seq1 = [1, 2, 3]
    seq2 = [1, 3]
    result = nw.align(seq1, seq2)
    assert result.length >= max(len(seq1), len(seq2))
    # Gaps are represented as 0.
    assert 0 in result.query or 0 in result.target


def test_blosum50_global():
    """Global alignment of two amino acid sequences."""
    scoring = Blosum50(gap_score=-8.0)
    nw = NeedlemanWunsch(scoring=scoring)
    seq1 = [13, 1, 15, 14, 6, 16, 3, 12, 4, 11]   # MAPFQSNKDL
    seq2 = [13, 11, 1, 15, 14, 7, 12, 17, 1, 1, 1, 2, 16, 10, 10]  # MLAPFEKTAAARSII
    result = nw.align(seq1, seq2)
    assert result.length >= max(len(seq1), len(seq2))
    assert isinstance(result.score, float)


def test_empty_sequences():
    scoring = SimpleMatch()
    nw = NeedlemanWunsch(scoring=scoring)
    result = nw.align([], [])
    assert result.length == 0
    assert result.score == 0.0
