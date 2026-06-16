"""Compute all-pairs alignment score matrices.

Given a collection of integer-encoded sequences and an aligner, produces a
symmetric distance/similarity matrix that can be fed into Neighbor-Joining
or other clustering algorithms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyseqalign.core.nw_affine import NeedlemanWunschAffine
from pyseqalign.scoring.protocols import ScoringFunction


@dataclass
class DistanceMatrix:
    """Symmetric distance/similarity matrix for a set of sequences.

    Attributes:
        labels: Ordered list of sequence labels (IDs).
        matrix: 2D numpy array where ``matrix[i, j]`` is the score between
            sequence *i* and sequence *j*.
        is_similarity: If ``True``, higher values mean more similar.
    """

    labels: list[str]
    matrix: np.ndarray
    is_similarity: bool = True

    @property
    def size(self) -> int:
        return len(self.labels)

    def get(self, label_a: str, label_b: str) -> float:
        """Look up the score between two sequences by label."""
        i = self.labels.index(label_a)
        j = self.labels.index(label_b)
        return float(self.matrix[i, j])

    def to_distance(self) -> DistanceMatrix:
        """Convert a similarity matrix to a distance matrix.

        Normalises similarities to [0, 1] then subtracts from 1.
        """
        if not self.is_similarity:
            return self
        mat = self.matrix.copy()
        diag_mask = np.eye(mat.shape[0], dtype=bool)
        off_diag = mat[~diag_mask]
        if off_diag.size == 0:
            return DistanceMatrix(
                labels=list(self.labels),
                matrix=np.zeros_like(mat),
                is_similarity=False,
            )
        min_val = off_diag.min()
        max_val = off_diag.max()
        rng = max_val - min_val
        if rng < 1e-12:
            dist = np.zeros_like(mat)
        else:
            dist = 1.0 - (mat - min_val) / rng
        np.fill_diagonal(dist, 0.0)
        return DistanceMatrix(
            labels=list(self.labels),
            matrix=dist,
            is_similarity=False,
        )

    def sub_matrix(self, label_subset: list[str]) -> DistanceMatrix:
        """Extract a sub-matrix for a subset of labels."""
        indices = [self.labels.index(lb) for lb in label_subset]
        sub = self.matrix[np.ix_(indices, indices)]
        return DistanceMatrix(
            labels=list(label_subset),
            matrix=sub.copy(),
            is_similarity=self.is_similarity,
        )


def compute_distance_matrix(
    sequences: dict[str, list[int]],
    scoring: ScoringFunction,
    gap_open: float = -2.5,
    gap_extend: float = -0.25,
) -> DistanceMatrix:
    """Compute all-pairs alignment scores for a set of sequences.

    Args:
        sequences: Mapping from sequence ID to integer-encoded sequence.
        scoring: Scoring function for the aligner.
        gap_open: Affine gap-open cost.
        gap_extend: Affine gap-extend cost.

    Returns:
        A ``DistanceMatrix`` with similarity scores (higher = more similar).
    """
    labels = sorted(sequences.keys())
    n = len(labels)
    mat = np.zeros((n, n), dtype=np.float64)

    aligner = NeedlemanWunschAffine(scoring, gap_open=gap_open, gap_extend=gap_extend)

    for i in range(n):
        seq_i = sequences[labels[i]]
        for j in range(i + 1, n):
            seq_j = sequences[labels[j]]
            result = aligner.align(seq_i, seq_j)
            mat[i, j] = result.score
            mat[j, i] = result.score
        # Self-alignment score on the diagonal.
        result_self = aligner.align(seq_i, seq_i)
        mat[i, i] = result_self.score

    return DistanceMatrix(labels=labels, matrix=mat, is_similarity=True)
