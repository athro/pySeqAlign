"""Consensus sequence construction from pairwise alignments.

When progressively aligning sequences via a guide tree, each merge step
aligns two groups and produces a *consensus* — a single representative
sequence that captures the most frequent element at each alignment position.
"""

from __future__ import annotations

from collections import Counter


def build_consensus(aligned_a: list[int], aligned_b: list[int]) -> list[int]:
    """Build a consensus from two aligned integer sequences.

    At each position, if both sequences have a non-gap element, the first
    (``aligned_a``) is preferred.  If one is a gap, the non-gap element is
    used.  If both are gaps, 0 (gap) is kept.

    Args:
        aligned_a: First aligned sequence (0 = gap).
        aligned_b: Second aligned sequence (0 = gap).

    Returns:
        Consensus sequence (same length as inputs).
    """
    assert len(aligned_a) == len(aligned_b), (
        f'Aligned sequences must have equal length, got {len(aligned_a)} and {len(aligned_b)}'
    )
    consensus: list[int] = []
    for a, b in zip(aligned_a, aligned_b):
        if a != 0:
            consensus.append(a)
        elif b != 0:
            consensus.append(b)
        else:
            consensus.append(0)
    return consensus


def profile_consensus(columns: list[list[int]]) -> list[int]:
    """Build a consensus from multiple aligned sequences (profile).

    For each alignment column, picks the most frequent non-gap element.
    Falls back to gap (0) if a column is all gaps.

    Args:
        columns: List of aligned sequences, all of the same length.
            Each is a ``list[int]`` where 0 represents a gap.

    Returns:
        A single consensus sequence.
    """
    if not columns:
        return []

    length = len(columns[0])
    for col in columns:
        assert len(col) == length, (
            f'All sequences must have equal length, got {len(col)} vs expected {length}'
        )

    consensus: list[int] = []
    for pos in range(length):
        elements = [col[pos] for col in columns if col[pos] != 0]
        if elements:
            # Most common non-gap element.
            counter = Counter(elements)
            consensus.append(counter.most_common(1)[0][0])
        else:
            consensus.append(0)
    return consensus
