"""Base types and protocol for ILP learning backends.

Defines the common interface that both Aleph and Popper backends implement.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass
class ILPTask:
    """Specification of an ILP learning task.

    An ILP task consists of background knowledge, positive/negative examples,
    and a language bias that constrains the hypothesis space.

    Attributes:
        background: Prolog clauses providing background knowledge (facts & rules).
        positive: Positive example facts (goals that should succeed).
        negative: Negative example facts (goals that should fail).
        bias: Language bias declarations (mode declarations for Aleph,
            head_pred/body_pred for Popper).
        settings: Additional ILP system settings as key-value pairs.
        work_dir: Directory for intermediate files.  Created automatically
            if not provided.
    """

    background: list[str] = field(default_factory=list)
    positive: list[str] = field(default_factory=list)
    negative: list[str] = field(default_factory=list)
    bias: list[str] = field(default_factory=list)
    settings: dict[str, str] = field(default_factory=dict)
    work_dir: Path | None = None


@dataclass
class LearnedProgram:
    """Result of an ILP learning run.

    Attributes:
        clauses: The learned Prolog clauses (rules).
        score: Quality score assigned by the learner (interpretation depends
            on the backend -- accuracy for Popper, coverage for Aleph).
        stats: Backend-specific statistics (runtime, nodes explored, etc.).
        raw_output: Full textual output from the ILP system.
    """

    clauses: list[str] = field(default_factory=list)
    score: float = 0.0
    stats: dict[str, object] = field(default_factory=dict)
    raw_output: str = ""

    @property
    def program_text(self) -> str:
        """Return the learned clauses as a single Prolog program string."""
        return "\n".join(self.clauses)


@runtime_checkable
class ILPLearner(Protocol):
    """Protocol that all ILP backends must implement."""

    def learn(self, task: ILPTask) -> LearnedProgram:
        """Run the ILP learner on the given task and return the result."""
        ...
