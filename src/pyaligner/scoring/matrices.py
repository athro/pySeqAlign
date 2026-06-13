"""Substitution matrices for amino acid sequence alignment.

Supports loading matrices dynamically from NCBI-format text files.
A set of commonly used BLOSUM and PAM matrices are bundled with
the package and can be loaded by name.

Example usage::

    # Load a bundled matrix by name
    scoring = SubstitutionMatrix.from_bundled("BLOSUM62")

    # Load from any NCBI-format file on disk
    scoring = SubstitutionMatrix.from_file("/path/to/my/MATRIX")

    # Download directly from NCBI FTP
    scoring = SubstitutionMatrix.from_ncbi("PAM120")

    # Legacy convenience alias (still works)
    scoring = Blosum50()
"""

from __future__ import annotations

from pathlib import Path
from typing import TextIO

# Standard one-letter amino acid codes, indexed 1..20 to match the legacy encoding.
# Index 0 is reserved for the gap character '-'.
AMINO_ACIDS = [
    "-",  # 0 -- gap
    "a",  # 1
    "r",  # 2
    "n",  # 3
    "d",  # 4
    "c",  # 5
    "q",  # 6
    "e",  # 7
    "g",  # 8
    "h",  # 9
    "i",  # 10
    "l",  # 11
    "k",  # 12
    "m",  # 13
    "f",  # 14
    "p",  # 15
    "s",  # 16
    "t",  # 17
    "w",  # 18
    "y",  # 19
    "v",  # 20
]

# Reverse lookup: one-letter code -> integer ID.
_AA_TO_ID: dict[str, int] = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}

# Directory containing bundled NCBI matrix files.
_MATRIX_DATA_DIR = Path(__file__).parent / "matrix_data"

# NCBI FTP base URL for substitution matrices.
_NCBI_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/blast/matrices"


def _parse_ncbi_matrix(source: TextIO) -> dict[tuple[int, int], float]:
    """Parse an NCBI-format substitution matrix into ``{(id_i, id_j): score}``.

    The format consists of:
    - Comment lines starting with ``#``
    - A header row of single-letter amino acid codes
    - Data rows: amino acid letter followed by integer scores

    Only the 20 standard amino acids (A R N D C Q E G H I L K M F P S T W Y V)
    are extracted; columns for ambiguity codes (B, Z, X) and stop (*) are ignored.
    """
    col_ids: list[int] = []
    matrix: dict[tuple[int, int], float] = {}

    for line in source:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        tokens = line.split()

        # Detect the header row: first token is a single letter that appears
        # in our amino acid alphabet (case-insensitive).
        if not col_ids:
            # Header row -- all tokens should be single letters.
            if all(len(t) == 1 for t in tokens):
                for t in tokens:
                    aa_id = _AA_TO_ID.get(t.lower(), -1)
                    col_ids.append(aa_id)
                continue
            # Some files start the header with a leading letter (the row label
            # coincides with the header).  Try treating the first token as a
            # letter and the rest as letters too.

        if not col_ids:
            continue

        # Data row: first token is the amino acid letter, rest are scores.
        row_aa = tokens[0].lower()
        row_id = _AA_TO_ID.get(row_aa, -1)
        if row_id <= 0:
            # Not a standard amino acid row (B, Z, X, *), skip.
            continue

        scores = tokens[1:]
        for col_idx, score_str in enumerate(scores):
            if col_idx >= len(col_ids):
                break
            col_id = col_ids[col_idx]
            if col_id <= 0:
                # Not a standard amino acid column, skip.
                continue
            val = float(score_str)
            matrix[(row_id, col_id)] = val
            matrix[(col_id, row_id)] = val

    return matrix


class SubstitutionMatrix:
    """A substitution matrix scoring function loaded from an NCBI-format file.

    Implements the ``ScoringFunction`` protocol so it can be passed directly
    to ``SmithWaterman`` or ``NeedlemanWunsch``.

    Args:
        matrix: Symmetric ``{(id_i, id_j): score}`` mapping.
        name: Human-readable name for the matrix (e.g. ``"BLOSUM62"``).
        gap_score: Score returned when either element is a gap (ID 0).
    """

    def __init__(
        self,
        matrix: dict[tuple[int, int], float],
        name: str = "custom",
        gap_score: float = -8.0,
    ) -> None:
        self._matrix = matrix
        self.name = name
        self._gap_score = gap_score

    def score(self, a: int, b: int) -> float:
        """Return the substitution score for element IDs *a* and *b*."""
        if a == 0 or b == 0:
            return self._gap_score
        return self._matrix.get((a, b), 0.0)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        gap_score: float = -8.0,
    ) -> SubstitutionMatrix:
        """Load a substitution matrix from an NCBI-format text file.

        Args:
            path: Path to the matrix file.
            gap_score: Score returned for gap characters.

        Example::

            scoring = SubstitutionMatrix.from_file("my_matrices/BLOSUM45")
        """
        path = Path(path)
        with open(path) as f:
            matrix = _parse_ncbi_matrix(f)
        return cls(matrix, name=path.stem, gap_score=gap_score)

    @classmethod
    def from_string(
        cls,
        text: str,
        name: str = "custom",
        gap_score: float = -8.0,
    ) -> SubstitutionMatrix:
        """Parse a substitution matrix from an NCBI-format string.

        Args:
            text: The matrix text in NCBI format.
            name: Name to assign to the matrix.
            gap_score: Score returned for gap characters.
        """
        import io

        matrix = _parse_ncbi_matrix(io.StringIO(text))
        return cls(matrix, name=name, gap_score=gap_score)

    @classmethod
    def from_bundled(
        cls,
        name: str,
        gap_score: float = -8.0,
    ) -> SubstitutionMatrix:
        """Load one of the bundled NCBI matrices by name.

        Available matrices: BLOSUM50, BLOSUM60, BLOSUM62, BLOSUM70, BLOSUM80,
        BLOSUM90, BLOSUM100, PAM50, PAM150, PAM200, PAM250.

        Args:
            name: Matrix name (case-insensitive), e.g. ``"BLOSUM62"`` or ``"pam250"``.
            gap_score: Score returned for gap characters.

        Raises:
            FileNotFoundError: If no bundled matrix with that name exists.

        Example::

            scoring = SubstitutionMatrix.from_bundled("PAM250")
        """
        path = _MATRIX_DATA_DIR / name.upper()
        if not path.exists():
            available = sorted(
                p.name for p in _MATRIX_DATA_DIR.iterdir()
                if p.is_file() and not p.name.startswith(".")
            )
            raise FileNotFoundError(
                f"No bundled matrix '{name}'. Available: {', '.join(available)}"
            )
        return cls.from_file(path, gap_score=gap_score)

    @classmethod
    def from_ncbi(
        cls,
        name: str,
        gap_score: float = -8.0,
    ) -> SubstitutionMatrix:
        """Download a substitution matrix directly from the NCBI FTP server.

        This fetches the matrix at runtime from
        ``https://ftp.ncbi.nlm.nih.gov/blast/matrices/<name>``.

        Args:
            name: Matrix name as it appears on the NCBI FTP server
                  (e.g. ``"BLOSUM45"``, ``"PAM120"``).
            gap_score: Score returned for gap characters.

        Raises:
            urllib.error.URLError: If the download fails.
        """
        import io
        import urllib.request

        url = f"{_NCBI_FTP_URL}/{name}"
        with urllib.request.urlopen(url) as resp:
            text = resp.read().decode("ascii")
        matrix = _parse_ncbi_matrix(io.StringIO(text))
        return cls(matrix, name=name, gap_score=gap_score)

    @classmethod
    def list_bundled(cls) -> list[str]:
        """Return the names of all bundled matrices."""
        if not _MATRIX_DATA_DIR.exists():
            return []
        return sorted(
            p.name for p in _MATRIX_DATA_DIR.iterdir()
            if p.is_file() and not p.name.startswith(".")
        )

    def __repr__(self) -> str:
        return f"SubstitutionMatrix(name={self.name!r}, gap_score={self._gap_score})"


# ---------------------------------------------------------------------------
# Legacy convenience aliases
# ---------------------------------------------------------------------------

# Keep the old hardcoded BLOSUM50 data for backward compatibility.
_BLOSUM50_RAW: dict[tuple[int, int], float] = {
    # a (1)
    (1, 1): 5, (1, 2): -2, (1, 3): -1, (1, 4): -2, (1, 5): -1,
    (1, 6): -1, (1, 7): -1, (1, 8): 0, (1, 9): -2, (1, 10): -1,
    (1, 11): -2, (1, 12): -1, (1, 13): -1, (1, 14): -3, (1, 15): -1,
    (1, 16): 1, (1, 17): 0, (1, 18): -3, (1, 19): -2, (1, 20): 0,
    # r (2)
    (2, 2): 7, (2, 3): -1, (2, 4): -2, (2, 5): -4,
    (2, 6): 1, (2, 7): 0, (2, 8): -3, (2, 9): 0, (2, 10): -4,
    (2, 11): -3, (2, 12): 3, (2, 13): -2, (2, 14): -3, (2, 15): -3,
    (2, 16): -1, (2, 17): -1, (2, 18): -3, (2, 19): -1, (2, 20): -3,
    # n (3)
    (3, 3): 7, (3, 4): 2, (3, 5): -2,
    (3, 6): 0, (3, 7): 0, (3, 8): 0, (3, 9): 1, (3, 10): -3,
    (3, 11): -4, (3, 12): 0, (3, 13): -2, (3, 14): -4, (3, 15): -2,
    (3, 16): 1, (3, 17): 0, (3, 18): -4, (3, 19): -2, (3, 20): -3,
    # d (4)
    (4, 4): 8, (4, 5): -4,
    (4, 6): 0, (4, 7): 2, (4, 8): -1, (4, 9): -1, (4, 10): -4,
    (4, 11): -4, (4, 12): -1, (4, 13): -4, (4, 14): -5, (4, 15): -1,
    (4, 16): 0, (4, 17): -1, (4, 18): -5, (4, 19): -3, (4, 20): -4,
    # c (5)
    (5, 5): 13, (5, 6): -3, (5, 7): -3, (5, 8): -3, (5, 9): -3, (5, 10): -2,
    (5, 11): -2, (5, 12): -3, (5, 13): -2, (5, 14): -2, (5, 15): -4,
    (5, 16): -1, (5, 17): -1, (5, 18): -5, (5, 19): -3, (5, 20): -1,
    # q (6)
    (6, 6): 7, (6, 7): 2, (6, 8): -2, (6, 9): 1, (6, 10): -3,
    (6, 11): -2, (6, 12): 2, (6, 13): 0, (6, 14): -4, (6, 15): -1,
    (6, 16): 0, (6, 17): -1, (6, 18): -1, (6, 19): -1, (6, 20): -3,
    # e (7)
    (7, 7): 6, (7, 8): -3, (7, 9): 0, (7, 10): -4,
    (7, 11): -3, (7, 12): 1, (7, 13): -2, (7, 14): -3, (7, 15): -1,
    (7, 16): -1, (7, 17): -1, (7, 18): -3, (7, 19): -2, (7, 20): -3,
    # g (8)
    (8, 8): 8, (8, 9): -2, (8, 10): -4,
    (8, 11): -4, (8, 12): -2, (8, 13): -3, (8, 14): -4, (8, 15): -2,
    (8, 16): 0, (8, 17): -2, (8, 18): -3, (8, 19): -3, (8, 20): -4,
    # h (9)
    (9, 9): 10, (9, 10): -4,
    (9, 11): -3, (9, 12): 0, (9, 13): -1, (9, 14): -1, (9, 15): -2,
    (9, 16): -1, (9, 17): -2, (9, 18): -3, (9, 19): 2, (9, 20): -4,
    # i (10)
    (10, 10): 5, (10, 11): 2, (10, 12): -3, (10, 13): 2, (10, 14): 0,
    (10, 15): -3, (10, 16): -3, (10, 17): -1, (10, 18): -3, (10, 19): -1, (10, 20): 4,
    # l (11)
    (11, 11): 5, (11, 12): -3, (11, 13): 3, (11, 14): 1,
    (11, 15): -4, (11, 16): -3, (11, 17): -1, (11, 18): -2, (11, 19): -1, (11, 20): 1,
    # k (12)
    (12, 12): 6, (12, 13): -2, (12, 14): -4,
    (12, 15): -1, (12, 16): 0, (12, 17): -1, (12, 18): -3, (12, 19): -2, (12, 20): -3,
    # m (13)
    (13, 13): 7, (13, 14): 0,
    (13, 15): -3, (13, 16): -2, (13, 17): -1, (13, 18): -1, (13, 19): 0, (13, 20): 1,
    # f (14)
    (14, 14): 8, (14, 15): -4, (14, 16): -3, (14, 17): -2,
    (14, 18): 1, (14, 19): 4, (14, 20): -1,
    # p (15)
    (15, 15): 10, (15, 16): -1, (15, 17): -1,
    (15, 18): -4, (15, 19): -3, (15, 20): -3,
    # s (16)
    (16, 16): 5, (16, 17): 2, (16, 18): -4, (16, 19): -2, (16, 20): -2,
    # t (17)
    (17, 17): 5, (17, 18): -3, (17, 19): -2, (17, 20): 0,
    # w (18)
    (18, 18): 15, (18, 19): 2, (18, 20): -3,
    # y (19)
    (19, 19): 8, (19, 20): -1,
    # v (20)
    (20, 20): 5,
}


class Blosum50:
    """BLOSUM50 substitution matrix scoring function (legacy convenience class).

    Element IDs follow the legacy encoding (1..20 for amino acids, 0 for gap).

    For new code, prefer ``SubstitutionMatrix.from_bundled("BLOSUM50")``.
    """

    def __init__(self, gap_score: float = -8.0) -> None:
        self._gap_score = gap_score
        # Build symmetric lookup.
        self._matrix: dict[tuple[int, int], float] = {}
        for (i, j), v in _BLOSUM50_RAW.items():
            self._matrix[(i, j)] = v
            self._matrix[(j, i)] = v

    def score(self, a: int, b: int) -> float:
        """Return BLOSUM50 score for element IDs *a* and *b*."""
        if a == 0 or b == 0:
            return self._gap_score
        return self._matrix.get((a, b), 0.0)
