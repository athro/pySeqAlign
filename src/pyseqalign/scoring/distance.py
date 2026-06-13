"""Distance-based scoring functions.

Includes the Nienhuys-Cheng atom distance from the legacy distances.pl, as well
as a simple identity-match scorer useful for testing.
"""

from __future__ import annotations


class SimpleMatch:
    """Simple identity-based scoring: +match_score for equal elements, +mismatch_score otherwise.

    Args:
        match_score: Score when elements are identical.
        mismatch_score: Score when elements differ.
        gap_score: Score for gap characters (element ID 0).
    """

    def __init__(
        self,
        match_score: float = 5.0,
        mismatch_score: float = -4.0,
        gap_score: float = -8.0,
    ) -> None:
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def score(self, a: int, b: int) -> float:
        if a == 0 or b == 0:
            return self.gap_score
        return self.match_score if a == b else self.mismatch_score


class AtomDistance:
    """Nienhuys-Cheng distance for structured atoms.

    This is a Python port of the recursive atom distance from the legacy
    distances.pl Prolog knowledge base.  It operates on structured
    representations where each atom is a tuple of ``(predicate, *args)`` and
    computes a normalised distance in [0, 1].

    For the integer-ID based interface used by the alignment algorithms, use
    ``AtomDistance`` with an *atom_store* mapping IDs to structured atoms.

    Args:
        atom_store: Mapping from integer element IDs to structured atoms
            (tuples).  ID 0 is reserved for gaps.
        gap_score: Score returned for gap characters.
        similarity: If ``True``, return ``1 - distance`` (similarity mode,
            matching the legacy ``sym`` mode).
    """

    def __init__(
        self,
        atom_store: dict[int, tuple] | None = None,
        gap_score: float = -1.0,
        similarity: bool = True,
    ) -> None:
        self.atom_store = atom_store or {}
        self.gap_score = gap_score
        self.similarity = similarity

    def score(self, a: int, b: int) -> float:
        """Return the (dis)similarity score between atom IDs *a* and *b*."""
        if a == 0 or b == 0:
            return self.gap_score

        atom_a = self.atom_store.get(a)
        atom_b = self.atom_store.get(b)

        if atom_a is None or atom_b is None:
            return self.gap_score

        dist = self._atom_distance(atom_a, atom_b)
        if self.similarity:
            return 1.0 - dist
        return dist

    def _atom_distance(self, a: tuple, b: tuple) -> float:
        """Recursive Nienhuys-Cheng distance between two structured atoms."""
        if a == b:
            return 0.0

        # Atoms must be tuples: (predicate, arg1, arg2, ...).
        if not isinstance(a, tuple) or not isinstance(b, tuple):
            return 1.0

        pred_a, *args_a = a
        pred_b, *args_b = b

        # Different predicate or arity => maximal distance.
        if pred_a != pred_b or len(args_a) != len(args_b):
            return 1.0

        if len(args_a) == 0:
            return 0.0

        total = sum(self._atom_distance(ai, bi) for ai, bi in zip(args_a, args_b))
        return total / (2 * len(args_a))
