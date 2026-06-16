"""The 'balloon' NLP running example from Karwath & Kersting (ILP 2006).

Five POS-tagged sentences (atoms ``pred(word)``, e.g. ``nn(balloon)``) are
aligned with the Nienhuys-Cheng atom distance and rendered as a relational
sequence logo. The two abstraction levels mirror the paper's "ground" (full
atom) and "abstract" (part-of-speech only) logos -- the abstract level makes
the shared syntactic skeleton (dt nn nn vbd prp in in dt ...) stand out even
where the words differ.

    python examples/balloon_logo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from pyseqalign.logo import relational_logo
from pyseqalign.msa import progressive_msa
from pyseqalign.scoring.distance import AtomDistance

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _data import Encoder, load_tsv  # noqa: E402

_HERE = Path(__file__).resolve().parent

LEVELS = [
    {'name': 'ground\n(pred(word))', 'symbol': lambda t: f'{t[0]}({t[1]})' if len(t) > 1 else t[0]},
    {'name': 'abstract\n(POS)', 'symbol': lambda t: t[0]},
]


def main() -> None:
    enc = Encoder()
    seqs = load_tsv(_HERE / 'data' / 'balloon.tsv', enc)
    scoring = AtomDistance(atom_store=enc.atom_store, gap_score=-0.5, similarity=True)
    # paper gap costs: opening 1.5, extension 0.5 (here as negative rewards, scaled).
    msa = progressive_msa(seqs, scoring, gap_open=-1.5, gap_extend=-0.5)
    rows = list(msa.aligned_sequences.values())
    out_dir = _HERE / 'output'; out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / 'logo_balloon.png'
    ic = relational_logo(rows, enc.atom_store, str(out),
                         title=f'Balloon NLP example ({len(seqs)} sentences)', levels=LEVELS)
    print(f'balloon: {len(seqs)} sentences, {len(rows[0]) if rows else 0} cols, '
          f'total IC {sum(ic):.1f} bits -> {out}')


if __name__ == '__main__':
    main()
