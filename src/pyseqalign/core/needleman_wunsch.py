"""Needleman-Wunsch global sequence alignment.

Translated from the legacy C implementation in pyAlign.c.
"""

from __future__ import annotations

from pyseqalign.core.alignment import AlignmentResult
from pyseqalign.core.smith_waterman import ScoringFunction


class NeedlemanWunsch:
    """Needleman-Wunsch global alignment.

    Args:
        scoring: A scoring function implementing the ``ScoringFunction`` protocol.
        gap_penalty: Cost applied when introducing a gap.  The scoring function is
            called with element ID ``0`` to represent a gap character.
    """

    def __init__(self, scoring: ScoringFunction, gap_penalty: float | None = None) -> None:
        self.scoring = scoring
        self._explicit_gap_penalty = gap_penalty

    @property
    def gap_penalty(self) -> float:
        """Return gap cost -- derived from scoring(0,0) when not set explicitly."""
        if self._explicit_gap_penalty is not None:
            return self._explicit_gap_penalty
        return self.scoring.score(0, 0)

    def align(self, seq1: list[int], seq2: list[int]) -> AlignmentResult:
        """Compute the optimal global alignment of *seq1* and *seq2*.

        Args:
            seq1: First input sequence (list of integer element IDs).
            seq2: Second input sequence.

        Returns:
            An ``AlignmentResult`` with aligned sequences, score, and length.
        """
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        gap = self.gap_penalty

        # Initialise F-matrix.
        f_matrix = [[0.0] * cols for _ in range(rows)]
        tb_matrix = [[-1.0] * cols for _ in range(rows)]

        # Fill border gaps.
        for i in range(1, rows):
            f_matrix[i][0] = gap * i
        for j in range(1, cols):
            f_matrix[0][j] = gap * j

        # Fill matrices.
        for i in range(1, rows):
            for j in range(1, cols):
                match = f_matrix[i - 1][j - 1] + self.scoring.score(seq1[i - 1], seq2[j - 1])
                delete = f_matrix[i - 1][j] + self.scoring.score(seq1[i - 1], 0)
                insert = f_matrix[i][j - 1] + self.scoring.score(0, seq2[j - 1])

                choices = [match, delete, insert]
                best = _argmax(choices)

                f_matrix[i][j] = choices[best]
                tb_matrix[i][j] = float(best)

        score = f_matrix[rows - 1][cols - 1]

        # Traceback.
        align1: list[int] = []
        align2: list[int] = []
        i = rows - 1
        j = cols - 1

        while i > 0 and j > 0:
            if tb_matrix[i][j] == 0.0:
                # Diagonal -- match/mismatch.
                i -= 1
                j -= 1
                align1.append(seq1[i])
                align2.append(seq2[j])
            elif tb_matrix[i][j] == 1.0:
                # Up -- gap in seq2.
                i -= 1
                align1.append(seq1[i])
                align2.append(0)
            else:
                # Left -- gap in seq1.
                j -= 1
                align1.append(0)
                align2.append(seq2[j])

        while i > 0:
            i -= 1
            align1.append(seq1[i])
            align2.append(0)

        while j > 0:
            j -= 1
            align1.append(0)
            align2.append(seq2[j])

        align1.reverse()
        align2.reverse()

        return AlignmentResult(
            query=align1,
            target=align2,
            score=score,
            length=len(align1),
        )


def _argmax(values: list[float]) -> int:
    best = 0
    for i in range(1, len(values)):
        if values[i] > values[best]:
            best = i
    return best
