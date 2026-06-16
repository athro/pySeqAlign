"""Protocols for scoring and gap cost functions."""

from __future__ import annotations

from typing import Protocol


class ScoringFunction(Protocol):
    """Protocol for scoring/distance functions used by alignment algorithms.

    Element ID 0 is reserved for the gap character.
    """

    def score(self, a: int, b: int) -> float:
        """Return the similarity score between elements *a* and *b*."""
        ...


class GapModel(Protocol):
    """Position-aware gap costs for affine/padding aligners."""

    def gap_open_cost(self, position: int, seq_len: int, element: int) -> float:
        """Return the gap-opening cost at the given position."""
        ...

    def gap_extend_cost(self, position: int, seq_len: int, element: int) -> float:
        """Return the gap-extension cost at the given position."""
        ...
