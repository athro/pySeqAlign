"""pyAligner -- Sequence alignment with Prolog-style distance functions and ILP learning."""

from pyaligner.core.smith_waterman import SmithWaterman
from pyaligner.core.needleman_wunsch import NeedlemanWunsch
from pyaligner.core.alignment import AlignmentResult, LocalAlignmentResult

__version__ = "0.1.0"

__all__ = [
    "SmithWaterman",
    "NeedlemanWunsch",
    "AlignmentResult",
    "LocalAlignmentResult",
]
