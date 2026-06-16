"""pySeqAlign -- Sequence alignment with Prolog-style distance functions and ILP learning."""

from pyseqalign.core.alignment import AlignmentResult, LocalAlignmentResult
from pyseqalign.core.needleman_wunsch import NeedlemanWunsch
from pyseqalign.core.smith_waterman import SmithWaterman

__version__ = "0.1.3"

__all__ = [
    "SmithWaterman",
    "NeedlemanWunsch",
    "AlignmentResult",
    "LocalAlignmentResult",
]
