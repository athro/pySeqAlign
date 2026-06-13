"""Example: Learning alignment rules with ILP.

Demonstrates how to use the AlignmentTaskBuilder to construct an ILP task
from labelled sequence pairs, and how to invoke either the Popper or Aleph
backend to learn Prolog rules.

This example builds a simple task where "similar" sequences share common
subsequences, then writes the ILP files to a temporary directory.  To
actually run the learner you need Popper or SWI-Prolog installed.
"""

import tempfile
from pathlib import Path

from pyaligner.learning import AlignmentTaskBuilder
from pyaligner.utils.helpers import encode_sequence

# ---- Build training data from amino acid sequences ----

# Positive pairs: sequences that *should* be considered similar.
positive_pairs = [
    ("HEAGAWGHEE", "PAWHEAE"),   # share AWG / AW motif
    ("MAPFQSNKDL", "MAPFEKTAA"), # share MAPF prefix
    ("GHEEHEAG",   "GHEEAWG"),   # share GHEE prefix
]

# Negative pairs: sequences that should *not* be considered similar.
negative_pairs = [
    ("HEAGAWGHEE", "MMMMMMMM"),
    ("MAPFQSNKDL", "WWWWWWWW"),
    ("GHEEHEAG",   "PPPPPPPP"),
]

builder = AlignmentTaskBuilder()

for s1, s2 in positive_pairs:
    builder.add_positive_pair(encode_sequence(s1), encode_sequence(s2))

for s1, s2 in negative_pairs:
    builder.add_negative_pair(encode_sequence(s1), encode_sequence(s2))

# Add amino acid facts as background knowledge.
for aa_id, aa_char in enumerate("arndcqeghilkmfpstwyv", start=1):
    builder.add_background_fact(f"amino_acid({aa_id}, '{aa_char}').")

# ---- Generate ILP files ----

# For Popper:
builder_popper = AlignmentTaskBuilder()
for s1, s2 in positive_pairs:
    builder_popper.add_positive_pair(encode_sequence(s1), encode_sequence(s2))
for s1, s2 in negative_pairs:
    builder_popper.add_negative_pair(encode_sequence(s1), encode_sequence(s2))
builder_popper.use_default_alignment_bias_popper()

with tempfile.TemporaryDirectory() as tmpdir:
    task = builder_popper.write_files(Path(tmpdir), name="similarity")
    print(f"Popper ILP files written to: {tmpdir}")
    print()

    for fname in ["bk.pl", "exs.pl", "bias.pl"]:
        fpath = Path(tmpdir) / fname
        print(f"--- {fname} ---")
        print(fpath.read_text()[:500])
        print()

# For Aleph:
builder_aleph = AlignmentTaskBuilder()
for s1, s2 in positive_pairs:
    builder_aleph.add_positive_pair(encode_sequence(s1), encode_sequence(s2))
for s1, s2 in negative_pairs:
    builder_aleph.add_negative_pair(encode_sequence(s1), encode_sequence(s2))
builder_aleph.use_default_alignment_bias_aleph()
builder_aleph.set_setting("i", "3")
builder_aleph.set_setting("noise", "0")

with tempfile.TemporaryDirectory() as tmpdir:
    task = builder_aleph.write_files(Path(tmpdir), name="similarity")
    print(f"Aleph ILP files written to: {tmpdir}")
    print()

    for fname in ["similarity.b", "similarity.f", "similarity.n"]:
        fpath = Path(tmpdir) / fname
        print(f"--- {fname} ---")
        print(fpath.read_text()[:500])
        print()

# ---- To actually run learning (requires Popper or SWI-Prolog) ----
#
# from pyaligner.learning.popper import PopperLearner
# learner = PopperLearner(timeout=60)
# result = learner.learn(task)
# print(result.program_text)
#
# Or with Aleph:
# from pyaligner.learning.aleph import AlephLearner
# learner = AlephLearner(induce_mode="induce")
# result = learner.learn(task)
# print(result.program_text)
