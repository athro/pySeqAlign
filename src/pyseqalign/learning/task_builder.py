"""Helpers for constructing ILP tasks from alignment data.

Converts alignment examples into the background knowledge, positive/negative
examples, and bias declarations needed by ILP learners.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from pyseqalign.learning.base import ILPTask


def _seq_to_prolog_list(seq: Sequence[int]) -> str:
    """Convert a Python integer list to a Prolog list string."""
    return "[" + ",".join(str(x) for x in seq) + "]"


class AlignmentTaskBuilder:
    """Build ILP tasks from alignment training data.

    This builder generates the Prolog facts and bias declarations required
    to learn distance/scoring predicates or alignment classification rules
    from labelled sequence pairs.

    Example usage::

        builder = AlignmentTaskBuilder()
        builder.add_positive_pair([1, 2, 3], [1, 2, 4], label="similar")
        builder.add_negative_pair([1, 2, 3], [10, 11, 12], label="similar")
        builder.add_background_fact("amino_acid(1, a).")
        task = builder.build()
    """

    def __init__(self) -> None:
        self._positive: list[str] = []
        self._negative: list[str] = []
        self._background: list[str] = []
        self._bias: list[str] = []
        self._settings: dict[str, str] = {}
        self._pair_id: int = 0

    # ------------------------------------------------------------------
    # Adding examples
    # ------------------------------------------------------------------

    def add_positive_pair(
        self,
        seq1: Sequence[int],
        seq2: Sequence[int],
        label: str = "similar",
    ) -> AlignmentTaskBuilder:
        """Add a positive example pair (sequences that *should* be related).

        Args:
            seq1: First sequence as integer IDs.
            seq2: Second sequence as integer IDs.
            label: Predicate name for the relation being learned.
        """
        self._pair_id += 1
        pid = self._pair_id
        s1 = _seq_to_prolog_list(seq1)
        s2 = _seq_to_prolog_list(seq2)
        self._positive.append(f"{label}(p{pid}).")
        self._background.append(f"seq1(p{pid},{s1}).")
        self._background.append(f"seq2(p{pid},{s2}).")
        return self

    def add_negative_pair(
        self,
        seq1: Sequence[int],
        seq2: Sequence[int],
        label: str = "similar",
    ) -> AlignmentTaskBuilder:
        """Add a negative example pair (sequences that should *not* be related)."""
        self._pair_id += 1
        pid = self._pair_id
        s1 = _seq_to_prolog_list(seq1)
        s2 = _seq_to_prolog_list(seq2)
        self._negative.append(f"{label}(p{pid}).")
        self._background.append(f"seq1(p{pid},{s1}).")
        self._background.append(f"seq2(p{pid},{s2}).")
        return self

    def add_background_fact(self, fact: str) -> AlignmentTaskBuilder:
        """Add a raw Prolog fact to the background knowledge."""
        self._background.append(fact)
        return self

    def add_background_rule(self, rule: str) -> AlignmentTaskBuilder:
        """Add a raw Prolog rule to the background knowledge."""
        self._background.append(rule)
        return self

    # ------------------------------------------------------------------
    # Bias / language declarations
    # ------------------------------------------------------------------

    def add_bias(self, declaration: str) -> AlignmentTaskBuilder:
        """Add a raw bias declaration (mode, determination, head_pred, etc.)."""
        self._bias.append(declaration)
        return self

    def set_setting(self, key: str, value: str) -> AlignmentTaskBuilder:
        """Set an ILP system parameter."""
        self._settings[key] = value
        return self

    def use_default_alignment_bias_aleph(
        self, label: str = "similar"
    ) -> AlignmentTaskBuilder:
        """Add default Aleph mode/determination declarations for sequence alignment.

        Sets up modes to learn rules about when two sequences are *label*.
        """
        self._bias.extend([
            f":- modeh(1, {label}(+pair)).",
            ":- modeb(*, seq1(+pair, -list)).",
            ":- modeb(*, seq2(+pair, -list)).",
            ":- modeb(*, member(-int, +list)).",
            ":- modeb(1, length(+list, -int)).",
            ":- modeb(1, sw_score(+pair, -float)).",
            ":- modeb(1, nw_score(+pair, -float)).",
            ":- modeb(1, score_above(+float, #float)).",
            f":- determination({label}/1, seq1/2).",
            f":- determination({label}/1, seq2/2).",
            f":- determination({label}/1, member/2).",
            f":- determination({label}/1, length/2).",
            f":- determination({label}/1, sw_score/2).",
            f":- determination({label}/1, nw_score/2).",
            f":- determination({label}/1, score_above/2).",
        ])
        return self

    def use_default_alignment_bias_popper(
        self, label: str = "similar"
    ) -> AlignmentTaskBuilder:
        """Add default Popper bias declarations for sequence alignment."""
        self._bias.extend([
            f"head_pred({label},1).",
            "body_pred(seq1,2).",
            "body_pred(seq2,2).",
            "body_pred(member,2).",
            "body_pred(length,2).",
            "body_pred(sw_score,2).",
            "body_pred(nw_score,2).",
            "body_pred(score_above,2).",
            "max_body(6).",
            "max_vars(6).",
        ])
        return self

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(self, work_dir: Path | None = None) -> ILPTask:
        """Construct the ILP task from the accumulated data."""
        return ILPTask(
            background=list(self._background),
            positive=list(self._positive),
            negative=list(self._negative),
            bias=list(self._bias),
            settings=dict(self._settings),
            work_dir=work_dir,
        )

    def write_files(self, directory: Path, name: str = "alignment") -> ILPTask:
        """Write .b, .f, .n, and bias files to *directory* and return the task.

        This produces files in the format expected by both Aleph and Popper.

        Args:
            directory: Target directory (created if it doesn't exist).
            name: Base filename (without extension).
        """
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        task = self.build(work_dir=directory)

        # Background knowledge
        bk_path = directory / f"{name}.b"
        bk_lines = list(task.bias) + [""] + list(task.background)
        if task.settings:
            for k, v in task.settings.items():
                bk_lines.insert(0, f":- set({k},{v}).")
        bk_path.write_text("\n".join(bk_lines) + "\n")

        # Positive examples
        (directory / f"{name}.f").write_text("\n".join(task.positive) + "\n")

        # Negative examples
        (directory / f"{name}.n").write_text("\n".join(task.negative) + "\n")

        # Popper-style separate bias file (also useful for Aleph via .b)
        (directory / "bias.pl").write_text("\n".join(task.bias) + "\n")

        # Popper exs.pl (combines pos/neg)
        popper_exs = []
        for p in task.positive:
            fact = p.rstrip(".")
            popper_exs.append(f"pos({fact}).")
        for n in task.negative:
            fact = n.rstrip(".")
            popper_exs.append(f"neg({fact}).")
        (directory / "exs.pl").write_text("\n".join(popper_exs) + "\n")

        # Popper bk.pl (background only, no bias)
        (directory / "bk.pl").write_text("\n".join(task.background) + "\n")

        return task
