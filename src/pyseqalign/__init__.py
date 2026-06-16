"""pySeqAlign -- Sequence alignment with Prolog-style distance functions and ILP learning."""

from pyseqalign.core.alignment import (
    AffineAlignmentResult,
    AlignmentResult,
    LocalAlignmentResult,
)
from pyseqalign.core.needleman_wunsch import NeedlemanWunsch
from pyseqalign.core.nw_affine import NeedlemanWunschAffine
from pyseqalign.core.smith_waterman import SmithWaterman

__version__ = "0.1.4"

__all__ = [
    "SmithWaterman",
    "NeedlemanWunsch",
    "NeedlemanWunschAffine",
    "AlignmentResult",
    "AffineAlignmentResult",
    "LocalAlignmentResult",
]
