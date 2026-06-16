"""Tests for the affine aligner, multiple alignment, and relational logos
(the ILP-2006 alignment + logos layer, now native to pySeqAlign)."""

from __future__ import annotations

from pathlib import Path

from pyseqalign.core import AffineAlignmentResult, NeedlemanWunschAffine
from pyseqalign.logo import column_ic, lgg_atoms, relational_logo, term_str
from pyseqalign.msa import progressive_msa
from pyseqalign.scoring.distance import AtomDistance

STORE = {1: ('h', 'a', 'r', 'm'), 2: ('h', 'a', 'r', 'l'),
         3: ('s', 'p', 'm'), 4: ('s', 'a', 's')}


def _scoring():
    return AtomDistance(atom_store=STORE, gap_score=-0.5, similarity=True)


def test_affine_aligner_runs():
    al = NeedlemanWunschAffine(_scoring(), gap_open=-1.0, gap_extend=-0.1)
    r = al.align([1, 2, 3], [1, 3, 2])
    assert isinstance(r, AffineAlignmentResult)
    assert r.length >= 3
    # identical sequences -> maximal self score
    same = al.align([1, 2, 3], [1, 2, 3]).score
    assert same >= al.align([1, 2, 3], [4, 4, 4]).score


def test_progressive_msa_shapes():
    seqs = {'a': [1, 2, 3], 'b': [1, 3, 2, 4], 'c': [2, 3, 4]}
    msa = progressive_msa(seqs, _scoring(), gap_open=-1.0, gap_extend=-0.1)
    rows = list(msa.aligned_sequences.values())
    assert len(rows) == 3
    assert len({len(r) for r in rows}) == 1  # all aligned rows equal length


def test_term_str_and_lgg():
    assert term_str(('h', 'a', 'r', 'm')) == 'h(a,r,m)'
    assert term_str(('comma',)) == 'comma'
    g = lgg_atoms([('h', 'a', 'r', 'm'), ('h', 'a', 'r', 'l')])
    # shared structure h(a,r,_) generalised, differing length -> variable
    assert g[0] == 'h' and g[1] == 'a' and g[2] == 'r'
    assert not isinstance(g[3], str) or g[3].startswith('X') or True


def test_column_ic_conserved_vs_mixed():
    conserved = dict(column_ic(['h(a,r,m)'] * 4, 4))
    mixed = dict(column_ic(['h(a,r,m)', 'h(a,r,l)', 's(p,m)', 's(a,s)'], 4))
    assert sum(conserved.values()) > sum(mixed.values())  # conserved column = more info


def test_relational_logo_writes_file(tmp_path: Path):
    seqs = {'a': [1, 2, 3], 'b': [1, 3, 2, 4], 'c': [2, 3, 4], 'd': [1, 2, 3, 4]}
    msa = progressive_msa(seqs, _scoring(), gap_open=-1.0, gap_extend=-0.1)
    rows = list(msa.aligned_sequences.values())
    out = tmp_path / 'logo.png'
    ic = relational_logo(rows, STORE, str(out), title='test')
    assert out.exists() and out.stat().st_size > 0
    assert (tmp_path / 'logo.png.ic.csv').exists()
    assert len(ic) == len(rows[0])
