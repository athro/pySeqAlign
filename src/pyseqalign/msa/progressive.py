"""Progressive multiple sequence alignment via a guide tree.

Follows the NJ guide tree merge order to progressively align sequences:
1. Compute all-pairs distance matrix.
2. Build NJ guide tree.
3. Follow merge order: align each pair, build consensus, feed into next step.

This mirrors the legacy ``multiAlign`` method in ``LogicAlign.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

from pyseqalign.core.nw_affine import NeedlemanWunschAffine
from pyseqalign.msa.consensus import build_consensus
from pyseqalign.msa.distance_matrix import DistanceMatrix, compute_distance_matrix
from pyseqalign.msa.guide_tree import TreeNode, neighbor_joining
from pyseqalign.scoring.protocols import ScoringFunction


@dataclass
class MSAResult:
    """Result of a multiple sequence alignment.

    Attributes:
        consensus: The final consensus sequence.
        aligned_sequences: Mapping from original sequence ID to its aligned form.
        merge_order: The order in which sequences were merged (for provenance).
        guide_tree: The NJ guide tree used.
    """

    consensus: list[int]
    aligned_sequences: dict[str, list[int]]
    merge_order: list[tuple[str, str, str]]
    guide_tree: TreeNode


def progressive_msa(
    sequences: dict[str, list[int]],
    scoring: ScoringFunction,
    gap_open: float = -2.5,
    gap_extend: float = -0.25,
    distance_matrix: DistanceMatrix | None = None,
) -> MSAResult:
    """Perform progressive multiple sequence alignment.

    Args:
        sequences: Mapping from sequence ID to integer-encoded sequence.
        scoring: Scoring function for the aligner.
        gap_open: Affine gap-open cost.
        gap_extend: Affine gap-extend cost.
        distance_matrix: Pre-computed distance matrix (optional; computed if
            not provided).

    Returns:
        An ``MSAResult`` with the final consensus and aligned sequences.
    """
    n = len(sequences)
    if n == 0:
        return MSAResult(
            consensus=[],
            aligned_sequences={},
            merge_order=[],
            guide_tree=TreeNode(name='empty'),
        )

    labels = sorted(sequences.keys())

    if n == 1:
        label = labels[0]
        seq = sequences[label]
        return MSAResult(
            consensus=list(seq),
            aligned_sequences={label: list(seq)},
            merge_order=[],
            guide_tree=TreeNode(name=label),
        )

    # Step 1: Compute distance matrix.
    if distance_matrix is None:
        distance_matrix = compute_distance_matrix(
            sequences, scoring, gap_open=gap_open, gap_extend=gap_extend
        )

    # Step 2: Build guide tree.
    guide_tree = neighbor_joining(distance_matrix)

    # Step 3: Get merge order.
    merges = guide_tree.merge_order()

    # Step 4: Progressive alignment.
    aligner = NeedlemanWunschAffine(scoring, gap_open=gap_open, gap_extend=gap_extend)

    # Pool of available sequences (original + consensus intermediates).
    pool: dict[str, list[int]] = {k: list(v) for k, v in sequences.items()}
    # Track which original sequences are represented by each pool entry.
    members: dict[str, list[str]] = {k: [k] for k in sequences}
    # Track the aligned form of each original sequence.
    aligned: dict[str, list[int]] = {}

    for merge_id, left_id, right_id in merges:
        seq_a = pool[left_id]
        seq_b = pool[right_id]

        result = aligner.align(seq_a, seq_b)
        consensus = build_consensus(result.query, result.target)

        # Propagate alignment to original members.
        # For the first merge of a leaf, its "aligned form" is the aligned output.
        # For later merges, the consensus from the previous step was already
        # gap-inserted, so the new alignment just adds more gaps.
        _propagate_gaps(aligned, members.get(left_id, []), seq_a, result.query)
        _propagate_gaps(aligned, members.get(right_id, []), seq_b, result.target)

        pool[merge_id] = consensus
        members[merge_id] = members.get(left_id, []) + members.get(right_id, [])

    # The final consensus.
    if merges:
        final_id = merges[-1][0]
        final_consensus = pool[final_id]
    else:
        final_consensus = list(next(iter(sequences.values())))

    # Make sure all aligned sequences have the same length.
    if aligned:
        max_len = max(len(v) for v in aligned.values())
        for k in aligned:
            if len(aligned[k]) < max_len:
                aligned[k].extend([0] * (max_len - len(aligned[k])))

    return MSAResult(
        consensus=final_consensus,
        aligned_sequences=aligned,
        merge_order=[(m, l, r) for m, l, r in merges],
        guide_tree=guide_tree,
    )


def _propagate_gaps(
    aligned: dict[str, list[int]],
    member_ids: list[str],
    old_seq: list[int],
    new_seq: list[int],
) -> None:
    """Update aligned forms of member sequences to match new gap insertions.

    When ``old_seq`` is re-aligned and becomes ``new_seq`` (with possibly more
    gaps), we need to propagate those gap insertions into all the original
    sequences that were already folded into ``old_seq``.
    """
    if not member_ids:
        return

    # Build gap-insertion map: how to transform old_seq → new_seq.
    # Walk both sequences; whenever new_seq has a gap that old_seq doesn't,
    # that's an insertion point.
    gap_map = _build_gap_map(old_seq, new_seq)

    for mid in member_ids:
        if mid not in aligned:
            # First time this sequence appears — its aligned form *is* new_seq
            # mapped through the original content.
            aligned[mid] = _apply_gap_map(old_seq, new_seq, gap_map, mid, aligned)
        else:
            # Already has an aligned form — insert new gaps at the right places.
            aligned[mid] = _insert_gaps(aligned[mid], gap_map)


def _build_gap_map(old_seq: list[int], new_seq: list[int]) -> list[int]:
    """Return a list of indices in new_seq where gaps were inserted.

    Walks old_seq and new_seq in parallel. Positions in new_seq that correspond
    to newly inserted gaps (not present in old_seq) are recorded.
    """
    insertions: list[int] = []
    oi = 0
    for ni in range(len(new_seq)):
        if oi < len(old_seq) and new_seq[ni] == old_seq[oi]:
            oi += 1
        else:
            # This position in new_seq is a gap insertion relative to old_seq.
            insertions.append(ni)
    return insertions


def _apply_gap_map(
    old_seq: list[int],
    new_seq: list[int],
    gap_positions: list[int],
    member_id: str,
    aligned: dict[str, list[int]],
) -> list[int]:
    """For a member seeing its first alignment, produce its aligned form.

    The member's original content is embedded in ``old_seq``. We map it to
    ``new_seq`` length by inserting gaps at the insertion positions.
    """
    # First-time: just use new_seq directly (the member *is* old_seq).
    return list(new_seq)


def _insert_gaps(seq: list[int], gap_positions: list[int]) -> list[int]:
    """Insert gap characters at specified positions into an existing sequence."""
    if not gap_positions:
        return seq
    result: list[int] = []
    gap_set = set(gap_positions)
    si = 0
    target_len = len(seq) + len(gap_positions)
    for pos in range(target_len):
        if pos in gap_set:
            result.append(0)
        else:
            if si < len(seq):
                result.append(seq[si])
                si += 1
            else:
                result.append(0)
    return result
