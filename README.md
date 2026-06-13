# pyAligner

Sequence alignment library with Prolog-style distance functions and ILP-based rule learning.

pyAligner provides Smith-Waterman (local) and Needleman-Wunsch (global) sequence alignment algorithms with pluggable distance/scoring functions. Distance functions can be defined natively in Python, via Prolog predicates, or *learned from data* using Inductive Logic Programming (ILP). This enables alignment of structured logical atoms beyond simple character or amino acid sequences.

## Features

- **Smith-Waterman** local alignment with k-best non-overlapping results
- **Needleman-Wunsch** global alignment
- **Prolog-based distance functions** via SWI-Prolog integration (optional)
- **Substitution matrices** -- BLOSUM (50, 60, 62, 70, 80, 90, 100) and PAM (50, 150, 200, 250) bundled; any NCBI-format matrix loadable from file or downloaded at runtime
- **Nienhuys-Cheng distance** for recursive structural comparison of logical atoms
- **ILP learning** of alignment rules and distance predicates:
  - **Aleph** backend -- classic ILP system (Srinivasan, 2001)
  - **Popper** backend -- modern ILP via learning from failures (Cropper & Morel, 2021)
- **Pure Python** core -- no C extension required (unlike the legacy version)

## Installation

### Prerequisites

- **Python >= 3.10** (tested on 3.10 -- 3.14)
- **SWI-Prolog** (optional, for Prolog-based scoring and ILP learning)
- **Clingo** (optional, for the Popper ILP backend)

```bash
# macOS (Homebrew)
brew install swi-prolog clingo

# Ubuntu / Debian
sudo apt install swi-prolog gringo
```

### Setting up with pyenv (recommended)

```bash
# Install a Python version
pyenv install 3.13.9

# Create a pyenv virtualenv
pyenv virtualenv 3.13.9 pyaligner-env

# Activate it for this project directory
cd pyAligner
pyenv local pyaligner-env

# Install all dependencies
pip install -r requirements.txt
pip install -e ".[learning,dev]"

# Register the Jupyter kernel so notebooks can use this environment
pip install ipykernel
python -m ipykernel install --user --name pyaligner-env --display-name "Python (pyAligner)"
```

### Setting up with venv (system Python)

```bash
cd pyAligner
python3 -m venv .venv
source .venv/bin/activate    # Linux / macOS
# .venv\Scripts\activate     # Windows

# Install all dependencies
pip install -r requirements.txt
pip install -e ".[learning,dev]"

# Register the Jupyter kernel
pip install ipykernel
python -m ipykernel install --user --name pyaligner-env --display-name "Python (pyAligner)"
```

### Optional dependencies

Install selectively if you don't need everything:

```bash
# Core only (no Prolog, no ILP)
pip install -e .

# Prolog-based distance functions (requires SWI-Prolog)
pip install -e ".[prolog]"

# Popper ILP backend (requires SWI-Prolog + Clingo)
pip install -e ".[popper]"

# All learning backends (Aleph + Popper)
pip install -e ".[learning]"

# Development tools (pytest, ruff, mypy)
pip install -e ".[dev]"
```

### Using requirements.txt

A `requirements.txt` is provided listing all Python dependencies (core,
Prolog, Popper, Jupyter, and development tools):

```bash
pip install -r requirements.txt
```

> **Note (Python 3.14):** Popper currently uses `pkg_resources`, which was
> removed from `setuptools>=81`. If you are on Python 3.14, install
> `setuptools<81` before Popper: `pip install 'setuptools<81'`.
>
> **Note:** The `popper-ilp` package on PyPI is outdated (v1.1.0, 2021).
> The project dependency points to the actively maintained
> [GitHub repository](https://github.com/logic-and-learning-lab/Popper) (v4.2.0+).

## Quick Start

```python
from pyaligner import SmithWaterman, NeedlemanWunsch
from pyaligner.scoring import Blosum50

# Local alignment with BLOSUM50 scoring
sw = SmithWaterman(scoring=Blosum50(), gap_penalty=-8.0)
results = sw.align(seq1, seq2, k=4, min_score=2.0)

for alignment in results:
    print(f"Score: {alignment.score}")
    print(f"  Query:  {alignment.query}")
    print(f"  Target: {alignment.target}")

# Global alignment
nw = NeedlemanWunsch(scoring=Blosum50(), gap_penalty=-8.0)
result = nw.align(seq1, seq2)
print(f"Score: {result.score}")
```

### Substitution Matrices

```python
from pyaligner.scoring import SubstitutionMatrix

# Load a bundled matrix by name
blosum62 = SubstitutionMatrix.from_bundled("BLOSUM62")
pam250   = SubstitutionMatrix.from_bundled("PAM250")

# Use with any aligner
sw = SmithWaterman(scoring=blosum62, gap_penalty=-8.0)

# Load from any NCBI-format file on disk
custom = SubstitutionMatrix.from_file("/path/to/my/MATRIX")

# Download directly from the NCBI FTP server at runtime
pam120 = SubstitutionMatrix.from_ncbi("PAM120")

# List all bundled matrices
print(SubstitutionMatrix.list_bundled())
# ['BLOSUM100', 'BLOSUM50', 'BLOSUM60', 'BLOSUM62', 'BLOSUM70',
#  'BLOSUM80', 'BLOSUM90', 'PAM150', 'PAM200', 'PAM250', 'PAM50']
```

### Using Prolog Distance Functions

```python
from pyaligner.prolog import PrologEngine

engine = PrologEngine()
engine.consult("my_distances.pl")

sw = SmithWaterman(scoring=engine, gap_penalty=-8.0)
results = sw.align(seq1, seq2)
```

### Learning Alignment Rules with ILP

```python
from pyaligner.learning import AlignmentTaskBuilder
from pyaligner.learning.popper import PopperLearner

# Build an ILP task from labelled sequence pairs.
builder = AlignmentTaskBuilder()
builder.add_positive_pair([1, 2, 3], [1, 2, 4], label="similar")
builder.add_positive_pair([5, 6, 7], [5, 6, 8], label="similar")
builder.add_negative_pair([1, 2, 3], [10, 11, 12], label="similar")
builder.use_default_alignment_bias_popper()

task = builder.build()

# Learn rules with Popper (modern, recommended).
learner = PopperLearner(timeout=60)
result = learner.learn(task)
print(result.program_text)
```

Or using the classic Aleph backend:

```python
from pyaligner.learning.aleph import AlephLearner

builder = AlignmentTaskBuilder()
# ... add examples ...
builder.use_default_alignment_bias_aleph()
task = builder.build()

learner = AlephLearner(induce_mode="induce")
result = learner.learn(task)
print(result.program_text)
```

## Project Structure

```
pyAligner/
├── src/
│   └── pyaligner/
│       ├── core/              # Alignment algorithms (SW, NW)
│       ├── prolog/            # Prolog engine integration
│       │   └── knowledge/     # Prolog knowledge bases (.pl files)
│       ├── scoring/           # Scoring/distance functions
│       │   └── matrix_data/  # Bundled NCBI substitution matrices
│       ├── learning/          # ILP learning backends
│       │   ├── base.py        # Common types & protocol
│       │   ├── task_builder.py # Convert alignment data -> ILP format
│       │   ├── aleph.py       # Aleph backend (classic ILP)
│       │   ├── popper.py      # Popper backend (modern ILP)
│       │   └── aleph_files/   # Bundled Aleph Prolog source
│       └── utils/             # Helpers and data structures
├── tests/
├── examples/
└── docs/
```

## ILP Backends

### Popper (recommended)

[Popper](https://github.com/logic-and-learning-lab/Popper) is a modern ILP system that learns from failures by combining ASP and Prolog. Advantages over Aleph:

- Learns optimal and recursive programs
- No metarules required
- Handles noise (MDL-based)
- Scales well (no grounding required)
- Native Python API

### Aleph (classic)

[Aleph](https://en.wikipedia.org/wiki/Aleph_(ILP)) is the classic ILP system by Ashwin Srinivasan, implementing Mode-Directed Inverse Entailment. The bundled `aleph_swi_ak.pl` is a SWI-Prolog compatible version ported from the legacy pySeqAlign codebase.

Supports multiple induction modes: `induce`, `induce_max`, `induce_cover`, `induce_tree`, `induce_features`, `induce_constraints`, `induce_incremental`, `induce_theory`.

### Other Modern ILP Systems

For reference, other notable systems in the field include:

- [PyILP](https://github.com/danyvarghese/PyILP) -- Python interface wrapping Aleph and Metagol for teaching
- [ILASP](https://github.com/ilaspltd/ILASP-releases) -- Answer Set Programming-based ILP
- [Metagol](https://github.com/metagol/metagol) -- Meta-Interpretive Learning
- [DeepStochLog](https://github.com/ML-KULeuven/deepstochlog) -- Neural-symbolic ILP combining logic and neural networks

## Background

This is a modern Python rewrite of the legacy pyAlign and pySeqAlign libraries. pyAlign used SWIG-wrapped C code with YAP Prolog bindings for alignment. pySeqAlign contained the Aleph ILP framework for learning rules from alignment examples. The new version is pure Python with optional SWI-Prolog integration via [Janus](https://www.swi-prolog.org/packs/list?p=janus) (the modern Python-Prolog bridge, replacing the older pyswip).

## License

MIT License. See [LICENSE](LICENSE) for details.
