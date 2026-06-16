#!/bin/bash
# Build the optional fast C++ affine-gap aligner (pybind11) for the active Python.
# Usage: PY=/path/to/python ./build_cpp_aligner.sh   (default: python3)
set -e
PY="${PY:-python3}"
HERE="$(cd "$(dirname "$0")" && pwd)"
PYINC=$("$PY" -c "import sysconfig;print(sysconfig.get_path('include'))")
PBINC=$("$PY" -c "import pybind11;print(pybind11.get_include())")
SUF=$("$PY" -c "import sysconfig;print(sysconfig.get_config_var('EXT_SUFFIX'))")
clang++ -O3 -std=c++14 -shared -undefined dynamic_lookup -fPIC \
  -I"$PYINC" -I"$PBINC" "$HERE/cpp_aligner.cpp" -o "$HERE/../cpp_aligner$SUF"
echo "built $HERE/../cpp_aligner$SUF"
