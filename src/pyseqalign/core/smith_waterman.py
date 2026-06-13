"""Smith-Waterman local sequence alignment.

Translated from the legacy C implementation in swAlign.c.
Computes the k best non-overlapping local alignments between two sequences.
"""

from __future__ import annotations

from typing import Protocol

from pyseqalign.core.alignment import KLocalAlignmentResults, LocalAlignmentResult


class ScoringFunction(Protocol):
    """Protocol for scoring/distance functions used by alignment algorithms."""

    def score(self, a: int, b: int) -> float:
        """Return the similarity score between elements *a* and *b*."""
        ...


class SmithWaterman:
    """Smith-Waterman local alignment.

    Args:
        scoring: A scoring function implementing the ``ScoringFunction`` protocol.
        gap_penalty: Cost applied when introducing a gap (should be positive;
            it is subtracted internally).
    """

    def __init__(self, scoring: ScoringFunction, gap_penalty: float = 8.0) -> None:
        self.scoring = scoring
        self.gap_penalty = gap_penalty

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def align(
        self,
        seq1: list[int],
        seq2: list[int],
        k: int = 1,
        cutoff: float = 0.0,
        min_score: float = 2.0,
    ) -> KLocalAlignmentResults:
        """Compute up to *k* best non-overlapping local alignments.

        Args:
            seq1: First input sequence (list of integer element IDs).
            seq2: Second input sequence.
            k: Maximum number of non-overlapping alignments to return.
            cutoff: Minimum cell value to keep in the F-matrix (default 0 for SW).
            min_score: Cells with score above this are considered trace start candidates.

        Returns:
            A ``KLocalAlignmentResults`` containing up to *k* alignments sorted
            by score descending.
        """
        if k == 0:
            return KLocalAlignmentResults()

        rows = len(seq1) + 1
        cols = len(seq2) + 1

        # Initialise F-matrix and traceback matrix.
        f_matrix = [[0.0] * cols for _ in range(rows)]
        traceback = [[(-10, -10)] * cols for _ in range(rows)]

        # Fill the matrices and collect high-scoring cells.
        max_traces: list[tuple[int, int]] = []
        d = self.gap_penalty

        for i in range(1, rows):
            for j in range(1, cols):
                match = f_matrix[i - 1][j - 1] + self.scoring.score(seq1[i - 1], seq2[j - 1])
                delete = f_matrix[i - 1][j] - d
                insert = f_matrix[i][j - 1] - d

                choices = [cutoff, match, delete, insert]
                best_idx = _argmax(choices)
                f_matrix[i][j] = choices[best_idx]

                if best_idx == 1:
                    traceback[i][j] = (i - 1, j - 1)
                elif best_idx == 2:
                    traceback[i][j] = (i - 1, j)
                elif best_idx == 3:
                    traceback[i][j] = (i, j - 1)

                if choices[best_idx] > min_score:
                    max_traces.append((i, j))

        # Generate all candidate traces (sorted by score descending).
        candidates = self._generate_traces(f_matrix, traceback, max_traces, rows, cols)

        # Select up to k non-overlapping alignments.
        selected: list[LocalAlignmentResult] = []
        for candidate in candidates:
            if len(selected) >= k:
                break
            if not any(self._overlaps(s, candidate) for s in selected):
                selected.append(candidate)

        return KLocalAlignmentResults(alignments=selected)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _generate_traces(
        f_matrix: list[list[float]],
        traceback: list[list[tuple[int, int]]],
        max_traces: list[tuple[int, int]],
        rows: int,
        cols: int,
    ) -> list[LocalAlignmentResult]:
        """Traceback from each high-scoring cell to produce alignment candidates."""
        results: list[LocalAlignmentResult] = []

        for end_i, end_j in max_traces:
            path_a: list[int] = []
            path_b: list[int] = []
            score = 0.0

            ci, cj = end_i, end_j
            while traceback[ci][cj] != (-10, -10):
                path_a.append(ci)
                path_b.append(cj)
                score += f_matrix[ci][cj]
                ci, cj = traceback[ci][cj]

            path_a.append(ci)
            path_b.append(cj)

            # Reverse to get start-to-end order.
            path_a.reverse()
            path_b.reverse()

            length = len(path_a)
            results.append(
                LocalAlignmentResult(
                    query_path=path_a,
                    target_path=path_b,
                    start_query=ci,
                    start_target=cj,
                    end_query=end_i,
                    end_target=end_j,
                    length=length,
                    score=score,
                )
            )

        # Sort by score descending.
        results.sort(key=lambda r: r.score, reverse=True)
        return results

    @staticmethod
    def _overlaps(a: LocalAlignmentResult, b: LocalAlignmentResult) -> bool:
        """Check whether two local alignments share any (i, j) cell."""
        cells_a = set(zip(a.query_path, a.target_path))
        cells_b = set(zip(b.query_path, b.target_path))
        return bool(cells_a & cells_b)


def _argmax(values: list[float]) -> int:
    """Return the index of the maximum value."""
    best = 0
    for i in range(1, len(values)):
        if values[i] > values[best]:
            best = i
    return best
