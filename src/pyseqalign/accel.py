"""Optional fast C++ affine-gap aligner (pybind11).

pySeqAlign's core is pure Python; this accelerated aligner is OPTIONAL. Build it
with ``src/pyseqalign/cpp/build_cpp_aligner.sh`` (needs a C++ compiler + pybind11).

It implements the same 3-matrix affine Needleman-Wunsch recurrence as a numeric
kernel over a dense score matrix, so results are interchangeable with a
pure-Python affine aligner. Exposes ``CppAligner(num_ids, gap_open, gap_extend)``
with ``set_matrix(flat_row_major)`` and ``align(q, t) -> AlignResult``
(``.score/.query/.target/.gap_opens/.gap_extensions/.length``).
"""
from __future__ import annotations


def cpp_available() -> bool:
    """True if the compiled extension is built and importable."""
    try:
        from . import cpp_aligner  # noqa: F401
        return True
    except Exception:
        return False


def load():
    """Return the compiled module (raises a helpful error if not built)."""
    try:
        from . import cpp_aligner
        return cpp_aligner
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            'pySeqAlign C++ aligner not built. Run '
            'src/pyseqalign/cpp/build_cpp_aligner.sh (needs a C++ compiler + pybind11).'
        ) from e
