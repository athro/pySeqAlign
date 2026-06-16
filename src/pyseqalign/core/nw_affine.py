"""Needleman-Wunsch global alignment with affine gap penalties.

Translated from the legacy C++ AlignerAffine::_align() implementation.
Uses three DP matrices (M, Ix, Iy) to distinguish gap-open from gap-extend.
"""

from __future__ import annotations

import numpy as np

from pyseqalign.core.alignment import AffineAlignmentResult
from pyseqalign.scoring.protocols import ScoringFunction

# Matrix indices.
_M = 0  # match/mismatch
_IX = 1  # gap in target (consuming query element)
_IY = 2  # gap in query (consuming target element)


class NeedlemanWunschAffine:
    """Needleman-Wunsch with affine gap penalties.

    Recurrences (similarity mode):
      M[i][j]  = score(q[i], t[j]) + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
      Ix[i][j] = max(M[i-1][j] + gap_open, Ix[i-1][j] + gap_extend, Iy[i-1][j] + gap_open)
      Iy[i][j] = max(M[i][j-1] + gap_open, Iy[i][j-1] + gap_extend, Ix[i][j-1] + gap_open)

    Args:
        scoring: Scoring function (element ID 0 = gap).
        gap_open: Cost for opening a new gap (should be negative for penalties).
        gap_extend: Cost for extending an existing gap (should be negative,
            typically less severe than gap_open).
    """

    def __init__(
        self,
        scoring: ScoringFunction,
        gap_open: float = -2.5,
        gap_extend: float = -0.25,
    ) -> None:
        self.scoring = scoring
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def align(self, seq1: list[int], seq2: list[int]) -> AffineAlignmentResult:
        """Compute the optimal global alignment with affine gap penalties.

        Args:
            seq1: Query sequence (list of integer element IDs).
            seq2: Target sequence.

        Returns:
            An ``AffineAlignmentResult`` with aligned sequences and gap statistics.
        """
        n = len(seq1)
        m = len(seq2)

        NEG_INF = -np.inf

        # F[k, i, j] for k in {M=0, Ix=1, Iy=2}
        F = np.full((3, n + 1, m + 1), NEG_INF, dtype=np.float64)
        # Traceback: B[k, i, j, :] = (from_k, from_i, from_j)
        B = np.full((3, n + 1, m + 1, 3), -1, dtype=np.int32)

        F[_M, 0, 0] = 0.0

        d = self.gap_open
        e = self.gap_extend

        # --- Border initialization: gaps along query (Ix column) ---
        for i0 in range(n):
            i = i0 + 1
            if i > 1:
                F[_IX, i, 0] = F[_IX, i - 1, 0] + e
            else:
                F[_IX, i, 0] = d
            B[_IX, i, 0] = [_IX, i - 1, 0]
            # M and Iy are -inf along this border (already set).

        # --- Border initialization: gaps along target (Iy row) ---
        for j0 in range(m):
            j = j0 + 1
            if j > 1:
                F[_IY, 0, j] = F[_IY, 0, j - 1] + e
            else:
                F[_IY, 0, j] = d
            B[_IY, 0, j] = [_IY, 0, j - 1]
            # M and Ix are -inf along this border (already set).

        # --- Main DP fill ---
        for i0 in range(n):
            i = i0 + 1
            for j0 in range(m):
                j = j0 + 1

                # Match/mismatch: diagonal transition.
                s = self.scoring.score(seq1[i - 1], seq2[j - 1])
                candidates_m = (
                    F[_M, i - 1, j - 1] + s,
                    F[_IX, i - 1, j - 1] + s,
                    F[_IY, i - 1, j - 1] + s,
                )
                best_k = _argmax3(candidates_m)
                F[_M, i, j] = candidates_m[best_k]
                B[_M, i, j] = [best_k, i - 1, j - 1]

                # Ix: gap in target (consume query[i], skip target).
                candidates_ix = (
                    F[_M, i - 1, j] + d,  # new gap
                    F[_IX, i - 1, j] + e,  # extend gap
                    F[_IY, i - 1, j] + d,  # new gap
                )
                best_k = _argmax3(candidates_ix)
                F[_IX, i, j] = candidates_ix[best_k]
                B[_IX, i, j] = [best_k, i - 1, j]

                # Iy: gap in query (skip query, consume target[j]).
                candidates_iy = (
                    F[_M, i, j - 1] + d,  # new gap
                    F[_IY, i, j - 1] + e,  # extend gap
                    F[_IX, i, j - 1] + d,  # new gap
                )
                best_k = _argmax3(candidates_iy)
                F[_IY, i, j] = candidates_iy[best_k]
                B[_IY, i, j] = [best_k, i, j - 1]

        # --- Find best endpoint ---
        end_scores = (F[_M, n, m], F[_IX, n, m], F[_IY, n, m])
        best_end = _argmax3(end_scores)
        score = end_scores[best_end]

        # --- Traceback ---
        align1, align2, gap_opens, gap_extensions = self._traceback(B, seq1, seq2, best_end, n, m)

        return AffineAlignmentResult(
            query=align1,
            target=align2,
            score=float(score),
            length=len(align1),
            gap_opens=gap_opens,
            gap_extensions=gap_extensions,
        )

    @staticmethod
    def _traceback(
        B: np.ndarray,
        seq1: list[int],
        seq2: list[int],
        start_k: int,
        start_i: int,
        start_j: int,
    ) -> tuple[list[int], list[int], int, int]:
        """Walk the traceback matrix to produce aligned sequences."""
        align1: list[int] = []
        align2: list[int] = []
        gap_opens = 0
        gap_extensions = 0

        k, i, j = start_k, start_i, start_j
        prev_k = -1

        while i > 0 or j > 0:
            from_k, from_i, from_j = int(B[k, i, j, 0]), int(B[k, i, j, 1]), int(B[k, i, j, 2])

            if from_i < 0 or from_j < 0:
                # Reached uninitialised border — shouldn't happen.
                break

            if k == _M:
                # Diagonal: match/mismatch.
                align1.append(seq1[i - 1])
                align2.append(seq2[j - 1])
            elif k == _IX:
                # Gap in target.
                align1.append(seq1[i - 1])
                align2.append(0)
                if prev_k != _IX:
                    gap_opens += 1
                else:
                    gap_extensions += 1
            else:  # _IY
                # Gap in query.
                align1.append(0)
                align2.append(seq2[j - 1])
                if prev_k != _IY:
                    gap_opens += 1
                else:
                    gap_extensions += 1

            prev_k = k
            k, i, j = from_k, from_i, from_j

        align1.reverse()
        align2.reverse()
        return align1, align2, gap_opens, gap_extensions


def _argmax3(vals: tuple[float, float, float]) -> int:
    """Return index of maximum among exactly three values."""
    if vals[0] >= vals[1]:
        return 0 if vals[0] >= vals[2] else 2
    return 1 if vals[1] >= vals[2] else 2
