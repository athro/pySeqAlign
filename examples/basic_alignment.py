"""Basic alignment example.

Reproduces the legacy swAlign_demo.py using the new pyseqalign API.
Aligns HEAGAWGHEE against PAWHEAE using Smith-Waterman with BLOSUM50.
"""

from pyseqalign import SmithWaterman, NeedlemanWunsch
from pyseqalign.scoring import Blosum50
from pyseqalign.utils.helpers import encode_sequence, decode_sequence

# Encode amino acid sequences to integer IDs.
seq1 = encode_sequence("HEAGAWGHEE")
seq2 = encode_sequence("PAWHEAE")

print("Sequences:")
print(f"  Seq1: HEAGAWGHEE  ->  {seq1}")
print(f"  Seq2: PAWHEAE     ->  {seq2}")
print()

# --- Smith-Waterman local alignment ---
scoring = Blosum50(gap_score=-8.0)
sw = SmithWaterman(scoring=scoring, gap_penalty=8.0)
results = sw.align(seq1, seq2, k=4, min_score=2.0)

print(f"Smith-Waterman: found {len(results)} alignment(s)")
for i, aln in enumerate(results):
    print(f"  Alignment {i + 1}:")
    print(f"    Score: {aln.score:.1f}")
    print(f"    Query region:  ({aln.start_query}, {aln.end_query})")
    print(f"    Target region: ({aln.start_target}, {aln.end_target})")
    print(f"    Length: {aln.length}")
print()

# --- Needleman-Wunsch global alignment ---
seq3 = encode_sequence("MAPFQSNKDL")
seq4 = encode_sequence("MLAPFEKTAAARSII")

print(f"Global alignment of MAPFQSNKDL vs MLAPFEKTAAARSII")
nw = NeedlemanWunsch(scoring=scoring)
result = nw.align(seq3, seq4)
print(f"  Score: {result.score:.1f}")
print(f"  Query:  {decode_sequence(result.query)}")
print(f"  Target: {decode_sequence(result.target)}")
print(f"  Length: {result.length}")
