"""Neighbor-Joining guide tree construction.

Translated from the legacy ``SimTree.py`` module. Produces a binary tree
whose leaf-traversal order defines the progressive alignment schedule.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyseqalign.msa.distance_matrix import DistanceMatrix

# ---------------------------------------------------------------------------
# Tree node
# ---------------------------------------------------------------------------


@dataclass
class TreeNode:
    """Node in a Neighbor-Joining guide tree.

    Attributes:
        name: Unique identifier (leaf label or internal merge ID).
        left: Left child (``None`` for leaves).
        right: Right child (``None`` for leaves).
        left_distance: Branch length to the left child.
        right_distance: Branch length to the right child.
    """

    name: str
    left: TreeNode | None = None
    right: TreeNode | None = None
    left_distance: float = 0.0
    right_distance: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return self.left is None and self.right is None

    def leaves(self) -> list[str]:
        """Return all leaf names reachable from this node (in-order)."""
        if self.is_leaf:
            return [self.name]
        result: list[str] = []
        if self.left is not None:
            result.extend(self.left.leaves())
        if self.right is not None:
            result.extend(self.right.leaves())
        return result

    def merge_order(self) -> list[tuple[str, str, str]]:
        """Return the progressive merge schedule as ``(new_id, left_id, right_id)`` tuples.

        Traverses the tree bottom-up so that every merge references only
        previously-available sequences or earlier merges.
        """
        merges: list[tuple[str, str, str]] = []
        self._collect_merges(merges)
        return merges

    def _collect_merges(self, acc: list[tuple[str, str, str]]) -> None:
        if self.is_leaf:
            return
        if self.left is not None:
            self.left._collect_merges(acc)
        if self.right is not None:
            self.right._collect_merges(acc)
        left_name = self.left.name if self.left else ''
        right_name = self.right.name if self.right else ''
        acc.append((self.name, left_name, right_name))


# ---------------------------------------------------------------------------
# Neighbor-Joining algorithm
# ---------------------------------------------------------------------------


def neighbor_joining(dist_matrix: DistanceMatrix) -> TreeNode:
    """Build a guide tree using the Neighbor-Joining algorithm.

    Expects a *distance* matrix (lower = more similar, diagonal = 0).

    Args:
        dist_matrix: A ``DistanceMatrix`` with ``is_similarity=False``.

    Returns:
        The root ``TreeNode`` of the NJ tree.
    """
    if dist_matrix.is_similarity:
        dist_matrix = dist_matrix.to_distance()

    labels = list(dist_matrix.labels)
    n = len(labels)

    if n == 0:
        return TreeNode(name='empty')
    if n == 1:
        return TreeNode(name=labels[0])
    if n == 2:
        d = dist_matrix.get(labels[0], labels[1])
        return TreeNode(
            name=f'_nj_{labels[0]}_{labels[1]}',
            left=TreeNode(name=labels[0]),
            right=TreeNode(name=labels[1]),
            left_distance=d / 2.0,
            right_distance=d / 2.0,
        )

    # Work with a mutable dense matrix.
    D = dist_matrix.matrix.copy()
    active = list(range(n))  # Indices into D that are still active.
    nodes: dict[int, TreeNode] = {i: TreeNode(name=labels[i]) for i in range(n)}

    next_id = n  # Counter for internal node indices.

    while len(active) > 2:
        k = len(active)

        # Compute r_i = sum of distances from i to all other active nodes.
        r = np.zeros(D.shape[0], dtype=np.float64)
        for idx in active:
            r[idx] = sum(D[idx, jdx] for jdx in active)

        # Find the pair (i, j) minimizing Q[i][j] = (k-2)*D[i][j] - r[i] - r[j].
        best_q = np.inf
        best_i, best_j = active[0], active[1]
        for ai in range(len(active)):
            for aj in range(ai + 1, len(active)):
                i, j = active[ai], active[aj]
                q = (k - 2) * D[i, j] - r[i] - r[j]
                if q < best_q:
                    best_q = q
                    best_i, best_j = i, j

        # Branch lengths from the new node to i and j.
        denom = 2.0 * (k - 2) if k > 2 else 2.0
        d_iu = D[best_i, best_j] / 2.0 + (r[best_i] - r[best_j]) / denom
        d_ju = D[best_i, best_j] - d_iu

        # Ensure non-negative branch lengths (NJ can produce small negatives).
        d_iu = max(d_iu, 0.0)
        d_ju = max(d_ju, 0.0)

        # Compute distances from the new node u to every remaining node m.
        # D[u][m] = (D[i][m] + D[j][m] - D[i][j]) / 2
        u = next_id
        next_id += 1

        # Expand D if needed.
        if u >= D.shape[0]:
            new_size = u + 1
            new_D = np.zeros((new_size, new_size), dtype=np.float64)
            new_D[: D.shape[0], : D.shape[1]] = D
            D = new_D

        for m in active:
            if m == best_i or m == best_j:
                continue
            d_um = (D[best_i, m] + D[best_j, m] - D[best_i, best_j]) / 2.0
            D[u, m] = d_um
            D[m, u] = d_um
        D[u, u] = 0.0

        # Create the new internal node.
        new_node = TreeNode(
            name=f'_nj_{u}',
            left=nodes[best_i],
            right=nodes[best_j],
            left_distance=d_iu,
            right_distance=d_ju,
        )
        nodes[u] = new_node

        # Update active set.
        active.remove(best_i)
        active.remove(best_j)
        active.append(u)

    # Join the last two remaining nodes.
    i, j = active[0], active[1]
    d = D[i, j]
    root = TreeNode(
        name='_nj_root',
        left=nodes[i],
        right=nodes[j],
        left_distance=d / 2.0,
        right_distance=d / 2.0,
    )
    return root
