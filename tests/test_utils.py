"""Tests for utility functions."""

from pyseqalign.utils.helpers import amino_acid_to_id, id_to_amino_acid, encode_sequence, decode_sequence


def test_amino_acid_roundtrip():
    for aa in "arndcqeghilkmfpstwyv":
        aid = amino_acid_to_id(aa)
        assert aid > 0
        assert id_to_amino_acid(aid) == aa


def test_encode_heagawghee():
    """Match the legacy example encoding."""
    result = encode_sequence("HEAGAWGHEE")
    assert result == [9, 7, 1, 8, 1, 18, 8, 9, 7, 7]


def test_decode_roundtrip():
    seq = "heagawghee"
    encoded = encode_sequence(seq)
    decoded = decode_sequence(encoded)
    assert decoded == seq


def test_gap_character():
    assert amino_acid_to_id("-") == 0
    assert id_to_amino_acid(0) == "-"


def test_unknown_returns_gap():
    assert amino_acid_to_id("X") == 0
