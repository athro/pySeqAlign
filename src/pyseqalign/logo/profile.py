"""Position-specific profiles for relational sequence logos.

A :class:`RelationalProfile` represents the frequency/probability of each
logical atom (identified by integer ID) at each alignment column.  This is
the relational analogue of a sequence logo profile used in bioinformatics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyseqalign.logo.probability import FreqDist, LidstoneProbDist, MLEProbDist


@dataclass(frozen=True)
class PositionProfile:
    """Profile for a single alignment column.

    Attributes:
        position: Column index (0-based).
        freqdist: Frequency distribution of atom IDs at this position.
        gap_count: Number of gap characters (ID 0) at this position.
        total_sequences: Total number of sequences contributing.
    """

    position: int
    freqdist: FreqDist
    gap_count: int
    total_sequences: int

    @property
    def gap_fraction(self) -> float:
        """Fraction of sequences that are gaps at this position."""
        if self.total_sequences == 0:
            return 0.0
        return self.gap_count / self.total_sequences

    @property
    def occupancy(self) -> float:
        """Fraction of sequences that are *not* gaps at this position."""
        return 1.0 - self.gap_fraction

    def prob_dist(
        self, smoothing: float = 0.0, bins: int | None = None
    ) -> MLEProbDist | LidstoneProbDist:
        """Build a probability distribution over non-gap atoms.

        Args:
            smoothing: Lidstone gamma. 0 gives MLE, >0 gives smoothed.
            bins: Total number of possible atom types (for smoothing).

        Returns:
            A probability distribution over atom IDs.
        """
        if smoothing <= 0.0:
            return MLEProbDist(self.freqdist)
        return LidstoneProbDist(self.freqdist, gamma=smoothing, bins=bins)

    def information_content(
        self,
        num_atom_types: int | None = None,
        smoothing: float = 0.0,
    ) -> float:
        """Compute information content (in bits) for this position.

        IC = log2(S) - H(pos)

        where ``S`` is the number of possible atom types and ``H`` is the
        Shannon entropy at this position.  Higher IC means more conservation.

        Args:
            num_atom_types: Total number of distinct atom types in the
                alphabet.  Defaults to the number of distinct atoms observed
                at this position.
            smoothing: Lidstone gamma for probability estimation.

        Returns:
            Information content in bits.
        """
        fd = self.freqdist
        if fd.total == 0:
            return 0.0

        if num_atom_types is None:
            num_atom_types = fd.num_bins

        if num_atom_types <= 1:
            return 0.0

        max_entropy = math.log2(num_atom_types)

        # Compute Shannon entropy.
        if smoothing > 0:
            pdist = LidstoneProbDist(fd, gamma=smoothing, bins=num_atom_types)
            entropy = 0.0
            for sample in fd.samples():
                p = pdist.prob(sample)
                if p > 0:
                    entropy -= p * math.log2(p)
            # Account for unseen atoms if bins > observed.
            unseen = num_atom_types - fd.num_bins
            if unseen > 0:
                p_unseen = pdist.prob(_UNSEEN_SENTINEL)
                if p_unseen > 0:
                    entropy -= unseen * p_unseen * math.log2(p_unseen)
        else:
            entropy = 0.0
            for sample in fd.samples():
                p = fd.freq(sample)
                if p > 0:
                    entropy -= p * math.log2(p)

        return max(0.0, max_entropy - entropy)


# Sentinel for querying unseen-sample probability in smoothed distributions.
_UNSEEN_SENTINEL = object()


class RelationalProfile:
    """Position-specific profile over an alignment of logical-atom sequences.

    Built from a set of aligned integer-ID sequences (where 0 = gap).
    Provides per-position frequency distributions and information content.
    """

    def __init__(
        self,
        aligned_sequences: list[list[int]],
        num_atom_types: int | None = None,
    ) -> None:
        """Create a profile from aligned sequences.

        Args:
            aligned_sequences: List of aligned sequences (equal length),
                using integer atom IDs with 0 = gap.
            num_atom_types: Total number of distinct atom types in the
                registry.  If ``None``, inferred from the data.
        """
        if not aligned_sequences:
            self._positions: list[PositionProfile] = []
            self._num_atom_types = num_atom_types or 0
            self._num_sequences = 0
            self._length = 0
            return

        self._num_sequences = len(aligned_sequences)
        self._length = len(aligned_sequences[0])

        for seq in aligned_sequences:
            if len(seq) != self._length:
                raise ValueError(
                    f'All sequences must have equal length, '
                    f'got {len(seq)} vs expected {self._length}'
                )

        # Build per-position frequency distributions.
        all_atoms: set[int] = set()
        self._positions = []
        for pos in range(self._length):
            fd = FreqDist()
            gap_count = 0
            for seq in aligned_sequences:
                atom_id = seq[pos]
                if atom_id == 0:
                    gap_count += 1
                else:
                    fd.inc(atom_id)
                    all_atoms.add(atom_id)
            self._positions.append(
                PositionProfile(
                    position=pos,
                    freqdist=fd,
                    gap_count=gap_count,
                    total_sequences=self._num_sequences,
                )
            )

        self._num_atom_types = num_atom_types if num_atom_types is not None else len(all_atoms)

    @property
    def length(self) -> int:
        """Number of alignment columns."""
        return self._length

    @property
    def num_sequences(self) -> int:
        """Number of sequences in the alignment."""
        return self._num_sequences

    @property
    def num_atom_types(self) -> int:
        """Total number of distinct atom types."""
        return self._num_atom_types

    def __len__(self) -> int:
        return self._length

    def __getitem__(self, index: int) -> PositionProfile:
        return self._positions[index]

    def __iter__(self):
        return iter(self._positions)

    def information_content(self, smoothing: float = 0.0) -> list[float]:
        """Compute information content at each position.

        Args:
            smoothing: Lidstone gamma for probability estimation.

        Returns:
            List of IC values (in bits), one per position.
        """
        return [
            pp.information_content(
                num_atom_types=self._num_atom_types,
                smoothing=smoothing,
            )
            for pp in self._positions
        ]

    def total_information(self, smoothing: float = 0.0) -> float:
        """Sum of information content across all positions."""
        return sum(self.information_content(smoothing=smoothing))

    def consensus(self) -> list[int]:
        """Return the consensus sequence (most frequent atom at each position).

        Gaps (0) are used where all sequences have a gap.
        """
        result: list[int] = []
        for pp in self._positions:
            if pp.freqdist.total == 0:
                result.append(0)
            else:
                result.append(pp.freqdist.max())
        return result

    def occupancy(self) -> list[float]:
        """Return per-position occupancy (1 - gap fraction)."""
        return [pp.occupancy for pp in self._positions]

    def __repr__(self) -> str:
        return (
            f'<RelationalProfile length={self._length} '
            f'sequences={self._num_sequences} '
            f'atom_types={self._num_atom_types}>'
        )
