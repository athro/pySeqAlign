"""Data structures for alignment results."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class AlignmentResult:
    """Result of a global (Needleman-Wunsch) alignment.

    Attributes:
        query: Aligned query sequence (with gaps represented as 0).
        target: Aligned target sequence (with gaps represented as 0).
        score: Alignment score.
        length: Length of the alignment.
    """

    query: list[int]
    target: list[int]
    score: float
    length: int


@dataclass
class AffineAlignmentResult(AlignmentResult):
    """Extended result from affine-gap alignment.

    Attributes:
        gap_opens: Number of gap-open events in both sequences combined.
        gap_extensions: Number of gap-extension events.
    """

    gap_opens: int = 0
    gap_extensions: int = 0


@dataclass
class LocalAlignmentResult:
    """Result of a single local (Smith-Waterman) alignment.

    Attributes:
        query_path: Indices along the query sequence in the alignment.
        target_path: Indices along the target sequence in the alignment.
        start_query: Start position in the query.
        start_target: Start position in the target.
        end_query: End position in the query.
        end_target: End position in the target.
        length: Length of the alignment path.
        score: Alignment score.
    """

    query_path: list[int]
    target_path: list[int]
    start_query: int
    start_target: int
    end_query: int
    end_target: int
    length: int
    score: float


@dataclass
class KLocalAlignmentResults:
    """Container for k non-overlapping local alignments.

    Attributes:
        alignments: List of local alignment results, sorted by score descending.
    """

    alignments: list[LocalAlignmentResult] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.alignments)

    def __getitem__(self, index: int) -> LocalAlignmentResult:
        return self.alignments[index]

    def __iter__(self):
        return iter(self.alignments)
