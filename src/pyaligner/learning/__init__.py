"""Inductive Logic Programming (ILP) backends for learning alignment rules.

This subpackage provides a common interface for learning scoring functions and
alignment rules from example alignments.  Two backends are supported:

- **Aleph** -- the classic ILP system (Srinivasan, 2001) via SWI-Prolog.
  Ported from the legacy pySeqAlign code.
- **Popper** -- a modern ILP system (Cropper & Morel, 2021) that learns from
  failures using ASP/SAT solvers.  Recommended for new projects.
"""

from pyaligner.learning.base import ILPTask, LearnedProgram, ILPLearner
from pyaligner.learning.task_builder import AlignmentTaskBuilder

__all__ = [
    "ILPTask",
    "LearnedProgram",
    "ILPLearner",
    "AlignmentTaskBuilder",
]
