"""Helper utilities for sequence encoding and decoding.

Ported from the legacy aminoAcids.pl mappings and C helper functions.
"""

from __future__ import annotations

from pyaligner.scoring.matrices import AMINO_ACIDS

# Build lookup tables from the canonical amino acid list.
_AA_TO_ID: dict[str, int] = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}
_ID_TO_AA: dict[int, str] = {idx: aa for idx, aa in enumerate(AMINO_ACIDS)}


def amino_acid_to_id(aa: str) -> int:
    """Convert a single-letter amino acid code to its integer ID.

    Returns 0 (gap) for unknown characters.
    """
    return _AA_TO_ID.get(aa.lower(), 0)


def id_to_amino_acid(element_id: int) -> str:
    """Convert an integer element ID back to its amino acid character."""
    return _ID_TO_AA.get(element_id, "-")


def encode_sequence(sequence: str) -> list[int]:
    """Encode an amino acid string into a list of integer IDs.

    Example::

        >>> encode_sequence("HEAGAWGHEE")
        [9, 7, 1, 8, 1, 18, 8, 9, 7, 7]
    """
    return [amino_acid_to_id(ch) for ch in sequence]


def decode_sequence(ids: list[int]) -> str:
    """Decode a list of integer IDs back into an amino acid string.

    Example::

        >>> decode_sequence([9, 7, 1, 8, 1, 18, 8, 9, 7, 7])
        'heagawghee'
    """
    return "".join(id_to_amino_acid(i) for i in ids)
