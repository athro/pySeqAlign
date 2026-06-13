"""Pairwise matching example.

Reproduces the legacy swAlign_demo2.py using the new pyseqalign API.
Generates random sequences and finds the best pairwise match for each.
"""

import random

from pyseqalign import SmithWaterman
from pyseqalign.scoring import Blosum50

random.seed(1)

# Generate 100 random amino acid sequences (IDs 1..20, lengths 5..9).
num_sequences = 100
sequences = [
    [random.randrange(1, 21) for _ in range(random.randrange(5, 10))]
    for _ in range(num_sequences)
]

print(f"Generated {len(sequences)} random sequences")
print()

scoring = Blosum50(gap_score=-8.0)
sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)

for i in range(len(sequences)):
    max_score = -100000.0
    best_j = -1
    best_result = None
    seq1 = sequences[i]

    for j in range(len(sequences)):
        if i == j:
            continue
        seq2 = sequences[j]
        results = sw.align(seq1, seq2, k=1, min_score=2.0)

        if results and results[0].score > max_score:
            max_score = results[0].score
            best_j = j
            best_result = results[0]

    if best_result:
        print(
            f"Seq {i:3d} best match = {best_j:3d}  "
            f"score = {max_score:7.1f}  "
            f"len(seq) = {len(seq1)}  "
            f"start({best_result.start_query},{best_result.start_target}) "
            f"end({best_result.end_query},{best_result.end_target})"
        )
