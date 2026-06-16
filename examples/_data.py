"""Tiny self-contained loaders for the bundled examples (no external deps).

Turns sequences of atom strings (e.g. ``h(a,r,m)``, ``nn(balloon)``) into the
int-id + ``atom_store`` representation the aligners and the logo renderer use.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

_ATOM = re.compile(r'^([a-z0-9_]+)(?:\((.*)\))?$', re.IGNORECASE)


def parse_atom(s: str) -> tuple:
    """``h(a,r,m)`` -> ``('h','a','r','m')``; ``comma`` -> ``('comma',)``.

    Flat atoms only (sufficient for the bundled SCOP/NLP examples)."""
    s = s.strip()
    m = _ATOM.match(s)
    if not m:
        return (s,)
    pred, args = m.group(1), m.group(2)
    if not args:
        return (pred,)
    return (pred, *[a.strip() for a in args.split(',')])


class Encoder:
    """Assigns stable int ids (>=1; 0 = gap) to distinct atoms; builds atom_store."""

    def __init__(self) -> None:
        self.ids: dict[tuple, int] = {}
        self.atom_store: dict[int, tuple] = {}

    def encode(self, atoms: list[str]) -> list[int]:
        out = []
        for a in atoms:
            t = parse_atom(a)
            if t not in self.ids:
                i = len(self.ids) + 1
                self.ids[t] = i
                self.atom_store[i] = t
            out.append(self.ids[t])
        return out


def load_scop_fold(xml_path: Path, enc: Encoder) -> dict[str, list[int]]:
    """Parse a secStruc fold_*.xml into {domain_id: [atom_ids]} using `enc`."""
    root = ET.parse(xml_path).getroot()
    seqs: dict[str, list[int]] = {}
    for seq in root.findall('sequence'):
        sid = seq.get('id')
        atoms = []
        for atom in seq.findall('atom'):
            sym = atom.find('symbol')
            text = (sym.text if sym is not None else atom.text) or ''
            if text.strip():
                atoms.append(text.strip())
        if sid and atoms:
            seqs[sid] = enc.encode(atoms)
    return seqs


def load_tsv(tsv_path: Path, enc: Encoder) -> dict[str, list[int]]:
    """Parse ``id <TAB> space-separated-atoms`` lines into {id: [atom_ids]}."""
    seqs: dict[str, list[int]] = {}
    for line in Path(tsv_path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        sid, atoms = line.split('\t')
        seqs[sid] = enc.encode(atoms.split())
    return seqs
