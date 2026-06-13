"""Popper ILP backend.

Popper (Cropper & Morel, 2021) is a modern ILP system that learns from
failures.  It combines answer set programming (ASP) with Prolog to search
the hypothesis space efficiently.

Popper is the **recommended** backend for new projects.  It can learn
recursive and optimal programs, handles noise, and does not require
metarules.

Install Popper with::

    pip install popper-ilp

Or from the repository::

    pip install git+https://github.com/logic-and-learning-lab/Popper@main

Requires SWI-Prolog (>=9.2.0) and Clingo (>=5.6.2).
"""

from __future__ import annotations

import importlib
import tempfile
from pathlib import Path

from pyseqalign.learning.base import ILPTask, LearnedProgram


class PopperLearner:
    """Popper ILP backend.

    Uses the Popper Python API to learn Prolog clauses from alignment
    examples.

    Args:
        timeout: Maximum learning time in seconds (default 300).
        max_literals: Maximum number of literals per clause (default 6).
        max_vars: Maximum number of variables per clause (default 6).
        noisy: Enable noise-tolerant learning (learns MDL programs).
        eval_timeout: Timeout for evaluating each hypothesis in seconds.
        extra_args: Additional keyword arguments passed to Popper's
            ``Settings`` constructor.
    """

    def __init__(
        self,
        timeout: int = 300,
        max_literals: int = 6,
        max_vars: int = 6,
        noisy: bool = False,
        eval_timeout: int = 10,
        **extra_args: object,
    ) -> None:
        self.timeout = timeout
        self.max_literals = max_literals
        self.max_vars = max_vars
        self.noisy = noisy
        self.eval_timeout = eval_timeout
        self.extra_args = extra_args
        self._check_available()

    @staticmethod
    def _check_available() -> None:
        """Verify that Popper is installed."""
        try:
            importlib.import_module("popper")
        except ImportError as exc:
            raise ImportError(
                "Popper is not installed. Install with:\n"
                "  pip install popper-ilp\n"
                "or:\n"
                "  pip install git+https://github.com/logic-and-learning-lab/Popper@main"
            ) from exc

    def learn(self, task: ILPTask) -> LearnedProgram:
        """Run Popper on the given task.

        Writes the task to Popper's expected directory format, calls the
        Popper learning loop, and returns the result.
        """
        from popper.loop import learn_solution
        from popper.util import Settings

        work_dir = task.work_dir or Path(tempfile.mkdtemp(prefix="pyseqalign_popper_"))
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        # Write Popper-format files.
        self._write_popper_files(task, work_dir)

        # Configure and run Popper.
        settings = Settings(
            kbpath=str(work_dir),
            timeout=self.timeout,
            eval_timeout=self.eval_timeout,
            **self.extra_args,
        )

        prog, score, stats = learn_solution(settings)

        # Format results.
        clauses: list[str] = []
        raw_output = ""
        if prog is not None:
            import io
            buf = io.StringIO()
            settings.print_prog_score(prog, score, file=buf)
            raw_output = buf.getvalue()
            clauses = [str(clause) for clause in prog]

        return LearnedProgram(
            clauses=clauses,
            score=score if score is not None else 0.0,
            stats=dict(stats) if stats else {},
            raw_output=raw_output,
        )

    # ------------------------------------------------------------------
    # File writing
    # ------------------------------------------------------------------

    @staticmethod
    def _write_popper_files(task: ILPTask, directory: Path) -> None:
        """Write Popper-format input files to *directory*."""
        # Background knowledge (bk.pl).
        (directory / "bk.pl").write_text("\n".join(task.background) + "\n")

        # Examples (exs.pl) -- Popper uses pos()/neg() wrappers.
        exs_lines: list[str] = []
        for p in task.positive:
            fact = p.rstrip(".")
            exs_lines.append(f"pos({fact}).")
        for n in task.negative:
            fact = n.rstrip(".")
            exs_lines.append(f"neg({fact}).")
        (directory / "exs.pl").write_text("\n".join(exs_lines) + "\n")

        # Bias (bias.pl) -- Popper's hypothesis language spec.
        (directory / "bias.pl").write_text("\n".join(task.bias) + "\n")


class PopperFallbackLearner:
    """Fallback Popper learner that calls Popper via subprocess.

    Use this when you cannot import Popper directly (e.g., version
    conflicts) but have the ``popper.py`` script available on PATH.

    Args:
        popper_cmd: Command to run Popper (default ``"python -m popper"``).
        timeout: Maximum learning time in seconds.
    """

    def __init__(
        self,
        popper_cmd: str = "python -m popper",
        timeout: int = 300,
    ) -> None:
        self.popper_cmd = popper_cmd
        self.timeout = timeout

    def learn(self, task: ILPTask) -> LearnedProgram:
        """Run Popper via subprocess."""
        import subprocess

        work_dir = task.work_dir or Path(tempfile.mkdtemp(prefix="pyseqalign_popper_"))
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        PopperLearner._write_popper_files(task, work_dir)

        try:
            result = subprocess.run(
                self.popper_cmd.split() + [str(work_dir)],
                capture_output=True,
                text=True,
                timeout=self.timeout,
            )
            raw_output = result.stdout + result.stderr
        except FileNotFoundError:
            raise RuntimeError(
                f"Popper not found via '{self.popper_cmd}'. "
                "Install with: pip install popper-ilp"
            )
        except subprocess.TimeoutExpired:
            return LearnedProgram(
                raw_output=f"Popper timed out after {self.timeout}s",
                stats={"timeout": True},
            )

        # Parse stdout for learned clauses.
        clauses = self._parse_output(raw_output)
        return LearnedProgram(
            clauses=clauses,
            raw_output=raw_output,
        )

    @staticmethod
    def _parse_output(raw_output: str) -> list[str]:
        """Extract clauses from Popper's stdout."""
        clauses: list[str] = []
        in_program = False
        for line in raw_output.splitlines():
            stripped = line.strip()
            if stripped.startswith("% ") or stripped == "":
                if in_program and stripped == "":
                    in_program = False
                continue
            if ":-" in stripped or (stripped.endswith(".") and not stripped.startswith("%")):
                # Looks like a Prolog clause.
                if not stripped.startswith("pos(") and not stripped.startswith("neg("):
                    clauses.append(stripped)
                    in_program = True
        return clauses
