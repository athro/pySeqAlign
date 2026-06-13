"""Aleph ILP backend.

Runs the Aleph ILP system (Srinivasan, 2001) via SWI-Prolog to learn
Prolog clauses from alignment examples.

Requires SWI-Prolog installed on the system and accessible via ``pyswip``
or the ``swipl`` command.

The bundled ``aleph_swi_ak.pl`` file is a SWI-Prolog compatible version
of Aleph 5, originally ported from YAP Prolog in the legacy pySeqAlign
codebase.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from pyaligner.learning.base import ILPTask, LearnedProgram


class AlephLearner:
    """Aleph ILP backend.

    Uses SWI-Prolog to run Aleph's ``induce/1`` on the provided task.

    Args:
        aleph_path: Path to the ``aleph_swi_ak.pl`` file.  Defaults to the
            bundled version shipped with pyaligner.
        swipl_cmd: Command to invoke SWI-Prolog (default ``"swipl"``).
        induce_mode: Aleph induction mode.  One of ``"induce"``,
            ``"induce_max"``, ``"induce_cover"``, ``"induce_tree"``,
            ``"induce_features"``, ``"induce_constraints"``,
            ``"induce_incremental"`` (default ``"induce"``).
        timeout: Maximum seconds for the SWI-Prolog process (default 300).
    """

    VALID_MODES = {
        "induce",
        "induce_max",
        "induce_cover",
        "induce_tree",
        "induce_features",
        "induce_constraints",
        "induce_incremental",
        "induce_theory",
    }

    def __init__(
        self,
        aleph_path: str | Path | None = None,
        swipl_cmd: str = "swipl",
        induce_mode: str = "induce",
        timeout: int = 300,
    ) -> None:
        if aleph_path is None:
            aleph_path = Path(__file__).parent / "aleph_files" / "aleph_swi_ak.pl"
        self.aleph_path = Path(aleph_path)
        self.swipl_cmd = swipl_cmd
        self.timeout = timeout

        if induce_mode not in self.VALID_MODES:
            raise ValueError(
                f"Unknown induce_mode '{induce_mode}'. "
                f"Valid modes: {sorted(self.VALID_MODES)}"
            )
        self.induce_mode = induce_mode

    def learn(self, task: ILPTask) -> LearnedProgram:
        """Run Aleph on the given task.

        Writes the task to temporary files, invokes SWI-Prolog with Aleph,
        and parses the output for learned clauses.
        """
        work_dir = task.work_dir or Path(tempfile.mkdtemp(prefix="pyaligner_aleph_"))
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        # Write Aleph-format files.
        bk_lines = []
        for k, v in task.settings.items():
            bk_lines.append(f":- set({k},{v}).")
        bk_lines.extend(task.bias)
        bk_lines.append("")
        bk_lines.extend(task.background)
        (work_dir / "task.b").write_text("\n".join(bk_lines) + "\n")
        (work_dir / "task.f").write_text("\n".join(task.positive) + "\n")
        (work_dir / "task.n").write_text("\n".join(task.negative) + "\n")

        # Construct SWI-Prolog script.
        aleph_abs = self.aleph_path.resolve()
        task_abs = (work_dir / "task").resolve()
        result_abs = (work_dir / "result.pl").resolve()

        script = (
            f":- consult('{aleph_abs}').\n"
            f":- read_all('{task_abs}').\n"
            f":- {self.induce_mode}.\n"
            f":- write_rules('{result_abs}').\n"
            f":- halt.\n"
        )
        script_path = work_dir / "run_aleph.pl"
        script_path.write_text(script)

        # Run SWI-Prolog.
        try:
            result = subprocess.run(
                [self.swipl_cmd, "-s", str(script_path)],
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(work_dir),
            )
            raw_output = result.stdout + result.stderr
        except FileNotFoundError:
            raise RuntimeError(
                f"SWI-Prolog not found at '{self.swipl_cmd}'. "
                "Install SWI-Prolog or set swipl_cmd to the correct path."
            )
        except subprocess.TimeoutExpired:
            return LearnedProgram(
                raw_output=f"Aleph timed out after {self.timeout}s",
                stats={"timeout": True},
            )

        # Parse results.
        clauses = self._parse_output(raw_output, result_abs)

        return LearnedProgram(
            clauses=clauses,
            score=self._extract_score(raw_output),
            stats=self._extract_stats(raw_output),
            raw_output=raw_output,
        )

    # ------------------------------------------------------------------
    # Parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_output(raw_output: str, result_file: Path) -> list[str]:
        """Extract learned clauses from Aleph output."""
        clauses: list[str] = []

        # Try reading the written rules file first.
        if result_file.exists():
            text = result_file.read_text()
            for line in text.strip().splitlines():
                line = line.strip()
                if line and not line.startswith("%"):
                    clauses.append(line)
            if clauses:
                return clauses

        # Fall back to parsing stdout for [Rule N] blocks.
        in_rule = False
        current: list[str] = []
        for line in raw_output.splitlines():
            if "[Rule" in line:
                in_rule = True
                current = []
                continue
            if in_rule:
                stripped = line.strip()
                if stripped == "":
                    if current:
                        clauses.append(" ".join(current))
                        current = []
                    in_rule = False
                else:
                    current.append(stripped)
        if current:
            clauses.append(" ".join(current))

        return clauses

    @staticmethod
    def _extract_score(raw_output: str) -> float:
        """Extract accuracy or coverage score from Aleph output."""
        for line in raw_output.splitlines():
            if "Accuracy" in line or "accuracy" in line:
                parts = line.split()
                for p in parts:
                    try:
                        return float(p.strip("()%,"))
                    except ValueError:
                        continue
        return 0.0

    @staticmethod
    def _extract_stats(raw_output: str) -> dict[str, object]:
        """Extract statistics from Aleph output."""
        stats: dict[str, object] = {}
        for line in raw_output.splitlines():
            if "clauses constructed" in line.lower():
                parts = line.split()
                for p in parts:
                    try:
                        stats["clauses_constructed"] = int(p)
                        break
                    except ValueError:
                        continue
            if "nodes explored" in line.lower() or "nodes visited" in line.lower():
                parts = line.split()
                for p in parts:
                    try:
                        stats["nodes_explored"] = int(p)
                        break
                    except ValueError:
                        continue
        return stats
