"""Best-effort build of the optional C++ affine-gap aligner.

pySeqAlign's core is pure Python. The C++ aligner (``pyseqalign.cpp_aligner``)
is an OPTIONAL accelerator (~100-270x faster, identical results). We try to
compile it at install time on any platform that has a C++ compiler + pybind11;
if that fails (no compiler, no pybind11, unsupported platform, ...) the install
STILL SUCCEEDS and the library transparently falls back to the pure-Python
aligner -- see ``pyseqalign.accel``.

This is why the project publishes an sdist (not a pure-Python wheel): pip builds
from source on the target machine, giving every install a chance to compile the
accelerator locally for its own Python/ABI/platform.
"""

from __future__ import annotations

import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Only build on POSIX (macOS/Linux); the hand-tuned flags below are GCC/Clang.
_CPP = 'src/pyseqalign/cpp/cpp_aligner.cpp'
ext_modules: list[Extension] = []
if sys.platform != 'win32':
    try:
        import pybind11
        ext_modules = [
            Extension(
                'pyseqalign.cpp_aligner',
                [_CPP],
                include_dirs=[pybind11.get_include()],
                language='c++',
                optional=True,  # setuptools won't fail the build if this ext won't compile
                extra_compile_args=['-O3', '-std=c++14'],
            )
        ]
    except Exception as exc:  # pybind11 missing -> skip the accelerator
        print(f'pyseqalign: skipping optional C++ aligner ({exc}); pure-Python fallback.')


class BestEffortBuildExt(build_ext):
    """Compile the accelerator if possible; never break the install if not."""

    def run(self) -> None:
        try:
            super().run()
        except Exception as exc:  # pragma: no cover - depends on build env
            print(f'pyseqalign: optional C++ aligner not built ({exc}); pure-Python fallback.')

    def build_extension(self, ext) -> None:
        try:
            super().build_extension(ext)
        except Exception as exc:  # pragma: no cover - depends on build env
            print(f'pyseqalign: optional C++ aligner not built ({exc}); pure-Python fallback.')


setup(ext_modules=ext_modules, cmdclass={'build_ext': BestEffortBuildExt})
