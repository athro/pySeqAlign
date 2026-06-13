"""Tests for scoring functions."""

import pytest

from pyseqalign.scoring import Blosum50, SimpleMatch, AtomDistance, SubstitutionMatrix


class TestBlosum50:
    def test_symmetric(self):
        b = Blosum50()
        assert b.score(1, 2) == b.score(2, 1)

    def test_diagonal_positive(self):
        b = Blosum50()
        for i in range(1, 21):
            assert b.score(i, i) > 0, f"Self-score for {i} should be positive"

    def test_gap_returns_gap_score(self):
        b = Blosum50(gap_score=-8.0)
        assert b.score(0, 5) == -8.0
        assert b.score(5, 0) == -8.0

    def test_known_values(self):
        b = Blosum50()
        # a-a should be 5
        assert b.score(1, 1) == 5.0
        # w-w should be 15
        assert b.score(18, 18) == 15.0
        # c-c should be 13
        assert b.score(5, 5) == 13.0


class TestSimpleMatch:
    def test_match(self):
        s = SimpleMatch(match_score=5.0, mismatch_score=-4.0)
        assert s.score(3, 3) == 5.0

    def test_mismatch(self):
        s = SimpleMatch(match_score=5.0, mismatch_score=-4.0)
        assert s.score(3, 7) == -4.0

    def test_gap(self):
        s = SimpleMatch(gap_score=-8.0)
        assert s.score(0, 5) == -8.0


class TestSubstitutionMatrix:
    """Tests for the dynamic SubstitutionMatrix loader."""

    def test_from_bundled_blosum62(self):
        m = SubstitutionMatrix.from_bundled("BLOSUM62")
        assert m.name == "BLOSUM62"
        # A-A = 4 in BLOSUM62
        assert m.score(1, 1) == 4.0
        # Symmetric
        assert m.score(1, 2) == m.score(2, 1)

    def test_from_bundled_case_insensitive(self):
        m = SubstitutionMatrix.from_bundled("blosum62")
        assert m.score(1, 1) == 4.0

    def test_from_bundled_pam250(self):
        m = SubstitutionMatrix.from_bundled("PAM250")
        # A-A = 2 in PAM250
        assert m.score(1, 1) == 2.0
        # W-W = 17 in PAM250
        assert m.score(18, 18) == 17.0
        # F-Y = 7 in PAM250
        assert m.score(14, 19) == 7.0

    def test_gap_score(self):
        m = SubstitutionMatrix.from_bundled("BLOSUM62", gap_score=-10.0)
        assert m.score(0, 5) == -10.0
        assert m.score(5, 0) == -10.0

    def test_all_bundled_matrices_load(self):
        for name in SubstitutionMatrix.list_bundled():
            if name == "__init__.py":
                continue
            m = SubstitutionMatrix.from_bundled(name)
            # Every matrix should have positive self-scores for standard AAs.
            for aa_id in range(1, 21):
                assert m.score(aa_id, aa_id) > 0, (
                    f"{name}: self-score for AA {aa_id} should be positive"
                )

    def test_all_bundled_symmetric(self):
        for name in SubstitutionMatrix.list_bundled():
            if name == "__init__.py":
                continue
            m = SubstitutionMatrix.from_bundled(name)
            for i in range(1, 21):
                for j in range(i + 1, 21):
                    assert m.score(i, j) == m.score(j, i), (
                        f"{name}: asymmetric score for ({i}, {j})"
                    )

    def test_from_string(self):
        text = """\
#  Test matrix
   A  R
A  4 -1
R -1  5
"""
        m = SubstitutionMatrix.from_string(text, name="test")
        assert m.name == "test"
        assert m.score(1, 1) == 4.0  # A-A
        assert m.score(2, 2) == 5.0  # R-R
        assert m.score(1, 2) == -1.0  # A-R
        assert m.score(2, 1) == -1.0  # R-A (symmetric)

    def test_from_bundled_not_found(self):
        with pytest.raises(FileNotFoundError, match="No bundled matrix"):
            SubstitutionMatrix.from_bundled("NONEXISTENT_MATRIX")

    def test_list_bundled(self):
        names = SubstitutionMatrix.list_bundled()
        assert "BLOSUM62" in names
        assert "PAM250" in names
        assert len(names) >= 11  # We bundled 11 matrices

    def test_repr(self):
        m = SubstitutionMatrix.from_bundled("BLOSUM80")
        assert "BLOSUM80" in repr(m)

    def test_blosum50_bundled_matches_legacy(self):
        """Verify the bundled BLOSUM50 file matches the hardcoded legacy data."""
        legacy = Blosum50()
        bundled = SubstitutionMatrix.from_bundled("BLOSUM50")
        for i in range(1, 21):
            for j in range(1, 21):
                assert legacy.score(i, j) == bundled.score(i, j), (
                    f"Mismatch at ({i}, {j}): "
                    f"legacy={legacy.score(i, j)}, bundled={bundled.score(i, j)}"
                )

    def test_use_with_smith_waterman(self):
        """SubstitutionMatrix works as a drop-in ScoringFunction."""
        from pyseqalign import SmithWaterman
        from pyseqalign.utils.helpers import encode_sequence

        scoring = SubstitutionMatrix.from_bundled("BLOSUM62")
        sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
        seq1 = encode_sequence("HEAGAWGHEE")
        seq2 = encode_sequence("PAWHEAE")
        results = sw.align(seq1, seq2, k=1, min_score=1.0)
        assert len(results) >= 1
        assert results[0].score > 0


class TestAtomDistance:
    def test_identical_atoms(self):
        store = {1: ("f", "a", "b"), 2: ("f", "a", "b")}
        ad = AtomDistance(atom_store=store, similarity=True)
        assert ad.score(1, 2) == 1.0

    def test_different_predicates(self):
        store = {1: ("f", "a"), 2: ("g", "a")}
        ad = AtomDistance(atom_store=store, similarity=True)
        assert ad.score(1, 2) == 0.0

    def test_partial_match(self):
        store = {1: ("f", "a", "b"), 2: ("f", "a", "c")}
        ad = AtomDistance(atom_store=store, similarity=True)
        s = ad.score(1, 2)
        assert 0.0 < s < 1.0

    def test_gap_handling(self):
        ad = AtomDistance(gap_score=-1.0)
        assert ad.score(0, 5) == -1.0
