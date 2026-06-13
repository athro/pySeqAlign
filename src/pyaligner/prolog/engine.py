"""SWI-Prolog engine wrapper for computing distance/scoring via Prolog predicates.

This replaces the legacy YAP Prolog bindings (yapBind.c) with SWI-Prolog
via the ``janus_swi`` package (the modern Python↔SWI-Prolog bridge).
The engine implements the ``ScoringFunction`` protocol so it can be passed
directly to ``SmithWaterman`` or ``NeedlemanWunsch``.

Requires the ``prolog`` optional dependency::

    pip install pyaligner[prolog]
"""

from __future__ import annotations

import importlib
from pathlib import Path


class PrologEngine:
    """SWI-Prolog based scoring function using Janus.

    Args:
        dist_mode: ``"sym"`` for similarity mode (1 - distance) or
            ``"dist"`` for raw distance.  Matches the legacy ``dist/7``
            predicate convention.
        distance_type: Name of the distance type predicate (default
            ``"atomDistance"``).
        distance_method: Distance method identifier (default ``"nc"`` for
            Nienhuys-Cheng).
        iteration: Iteration parameter passed to ``dist/7`` (default ``0``).
        gap_default: Default score returned for gap characters.
    """

    def __init__(
        self,
        dist_mode: str = "sym",
        distance_type: str = "atomDistance",
        distance_method: str = "nc",
        iteration: int = 0,
        gap_default: float = -8.0,
    ) -> None:
        self.dist_mode = dist_mode
        self.distance_type = distance_type
        self.distance_method = distance_method
        self.iteration = iteration
        self.gap_default = gap_default

        # Lazy import so the package works without janus_swi installed.
        try:
            self._janus = importlib.import_module("janus_swi")
        except ImportError as exc:
            raise ImportError(
                "janus_swi is required for Prolog support. "
                "Install with: pip install pyaligner[prolog]"
            ) from exc

    def consult(self, path: str | Path) -> None:
        """Load a Prolog file into the engine."""
        self._janus.consult(str(path))

    def assert_fact(self, term: str) -> None:
        """Assert a Prolog fact."""
        self._janus.query_once(f"assert({term})")

    def retract_fact(self, term: str) -> None:
        """Retract a Prolog fact."""
        self._janus.query_once(f"retract({term})")

    def call(self, goal: str) -> list[dict]:
        """Execute a Prolog goal and return all solution bindings."""
        return list(self._janus.query(goal))

    def call_once(self, goal: str) -> dict:
        """Execute a Prolog goal and return the first solution."""
        return self._janus.query_once(goal)

    def score(self, a: int, b: int) -> float:
        """Compute the score between element IDs *a* and *b* via Prolog.

        Calls the ``dist/7`` predicate:
        ``dist(Mode, Type, Method, Iteration, A, B, D)``
        and returns the bound value of ``D``.
        """
        if a == 0 or b == 0:
            return self.gap_default

        query = (
            f"dist({self.dist_mode},{self.distance_type},"
            f"{self.distance_method},{self.iteration},{a},{b},D)"
        )
        result = self._janus.query_once(query)
        if result["truth"]:
            return float(result["D"])
        return -100.0

    def consult_knowledge_base(self) -> None:
        """Load the bundled Prolog knowledge base files (amino acids, BLOSUM50, distances)."""
        kb_dir = Path(__file__).parent / "knowledge"
        for pl_file in ["amino_acids.pl", "distances.pl", "blosum50.pl"]:
            path = kb_dir / pl_file
            if path.exists():
                self.consult(path)
