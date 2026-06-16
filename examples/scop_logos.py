"""Reproduce the SCOP relational sequence logos from Karwath & Kersting (ILP 2006).

For each SCOP fold (secondary-structure-element sequences), build a multiple
alignment with the Nienhuys-Cheng atom distance and render the relational logo
(stacked information-content glyphs + per-column least-general generalisation),
as in Figures 3-6 of the paper. Pure pySeqAlign -- no learning involved.

    python examples/scop_logos.py            # all bundled folds -> examples/output/
    python examples/scop_logos.py fold_003_001

SCOP data (c) SCOP/ASTRAL; see examples/data/README.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

from pyseqalign.logo import relational_logo
from pyseqalign.msa import progressive_msa
from pyseqalign.scoring.distance import AtomDistance

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _data import Encoder, load_scop_fold  # noqa: E402

_HERE = Path(__file__).resolve().parent
FOLD_NAMES = {
    'fold_003_001': 'c.1 TIM beta/alpha-barrel',
    'fold_003_002': 'c.2 NAD(P)-binding Rossmann fold',
    'fold_002_001': 'b.1 Immunoglobulin-like beta-sandwich',
    'fold_002_032': 'b.32 SH3-like barrel',
    'fold_001_002': 'a.2 Long alpha-hairpin',
}

# Two abstraction levels: the full SSE atom, and just the element type (h vs s).
LEVELS = [
    {'name': 'element', 'symbol': lambda t: f'{t[0]}({",".join(t[1:])})' if len(t) > 1 else t[0]},
    {'name': 'type', 'symbol': lambda t: t[0]},
]


def run_fold(xml: Path, out_dir: Path) -> None:
    enc = Encoder()
    seqs = load_scop_fold(xml, enc)
    if len(seqs) < 2:
        print(f'skip {xml.name}: <2 sequences'); return
    scoring = AtomDistance(atom_store=enc.atom_store, gap_score=-0.5, similarity=True)
    msa = progressive_msa(seqs, scoring, gap_open=-1.0, gap_extend=-0.1)
    rows = list(msa.aligned_sequences.values())
    title = f'{FOLD_NAMES.get(xml.stem, xml.stem)}  ({len(seqs)} domains)'
    out = out_dir / f'logo_{xml.stem}.png'
    ic = relational_logo(rows, enc.atom_store, str(out), title=title, levels=LEVELS)
    print(f'{xml.stem}: {len(seqs)} domains, {len(rows[0]) if rows else 0} cols, '
          f'total IC {sum(ic):.1f} bits -> {out.name}')


def main() -> None:
    out_dir = _HERE / 'output'
    out_dir.mkdir(parents=True, exist_ok=True)
    data = _HERE / 'data' / 'scop'
    wanted = sys.argv[1:]
    folds = ([data / f'{w}.xml' for w in wanted] if wanted
             else sorted(data.glob('fold_*.xml')))
    for xml in folds:
        run_fold(xml, out_dir)
    print(f'\nlogos -> {out_dir}/')


if __name__ == '__main__':
    main()
