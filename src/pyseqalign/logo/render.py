"""Relational sequence logos (Karwath & Kersting, ILP 2006).

Renders the information-content logo of a multiple alignment of *structured*
atoms. Atoms are represented as nested tuples ``(predicate, arg1, arg2, ...)``
-- the same representation used by :class:`pyseqalign.scoring.distance.AtomDistance`
-- so this module has no dependency on the boosting/learning layer.

Per-column information content follows the paper's Gorodkin/KL form
    I_i = sum_k q_ik * log2(q_ik / p_k)         (gap prior p_- = 1)
with per-symbol contributions I_ik = q_ik*log2(q_ik/p_k) giving the stacked
glyph heights; symbols rarer than expected (I_ik < 0) are drawn upside-down.

Typical use::

    from pyseqalign.logo.render import relational_logo
    relational_logo(aligned_rows, atom_store, 'logo.png', title='Fold c.1')

where ``aligned_rows`` is a list of equal-or-ragged int-id rows (gap = 0) from
a multiple alignment, and ``atom_store`` maps id -> structured tuple.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable

_PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
            '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


# --------------------------------------------------------------------------- #
# structured-atom helpers (tuples; gap = None)
# --------------------------------------------------------------------------- #
def term_str(t) -> str:
    """Render a tuple atom ``('h','a','r','m')`` as ``h(a,r,m)`` (nullary: ``h``)."""
    if t is None:
        return '-'
    if isinstance(t, tuple):
        pred, *args = t
        if not args:
            return str(pred)
        return f'{pred}({",".join(term_str(a) for a in args)})'
    return str(t)


def _freeze(t):
    return (t[0], tuple(_freeze(x) for x in t[1:])) if isinstance(t, tuple) else t


class _Var:
    __slots__ = ('name',)

    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:
        return self.name


def _lgg_pair(a, b, table: dict, counter: list):
    """Anti-unify two (possibly nested) tuple atoms; consistent vars via table."""
    if isinstance(a, tuple) and isinstance(b, tuple) and a and b \
            and a[0] == b[0] and len(a) == len(b):
        return (a[0], *[_lgg_pair(x, y, table, counter) for x, y in zip(a[1:], b[1:])])
    if a == b and not isinstance(a, tuple):
        return a
    key = (_freeze(a), _freeze(b))
    if key not in table:
        counter[0] += 1
        table[key] = _Var(f'X{counter[0]}')
    return table[key]


def lgg_atoms(atoms: list) -> object | None:
    """Least-general generalisation (anti-unification) of several tuple atoms."""
    atoms = [a for a in atoms if a is not None]
    if not atoms:
        return None
    acc = atoms[0]
    for a in atoms[1:]:
        acc = _lgg_pair(acc, a, {}, [0])
    return acc


# --------------------------------------------------------------------------- #
# information content
# --------------------------------------------------------------------------- #
def column_ic(symbols: list[str], alphabet_size: int) -> list[tuple[str, float]]:
    """Per-symbol IC contribution I_ik = q_ik*log2(q_ik/p_k) for one column.

    p_k = 1/alphabet_size for non-gap symbols; gap '-' excluded (prior p_- = 1).
    """
    n = len(symbols)
    if n == 0:
        return []
    counts = Counter(symbols)
    p = 1.0 / max(1, alphabet_size)
    out = []
    for sym, c in counts.items():
        if sym == '-':
            continue
        q = c / n
        out.append((sym, q * math.log2(q / p)))
    return out


# --------------------------------------------------------------------------- #
# drawing primitives
# --------------------------------------------------------------------------- #
def _glyph(ax, x, y0, sym, h, color, flip):
    if h <= 1e-6:
        return
    fs = max(3, min(22, h * 22))
    ax.text(x, y0 + h / 2, sym, fontsize=fs, ha='center', va='center',
            color=color, family='monospace', rotation=180 if flip else 0, clip_on=True)


def draw_stack_logo(ax, per_col, title, ymax):
    """per_col[i] = [(symbol, signed_height, color)] -- stacked IC logo for column i."""
    for x, stack in enumerate(per_col):
        up = sorted([s for s in stack if s[1] >= 0], key=lambda s: s[1])
        dn = sorted([s for s in stack if s[1] < 0], key=lambda s: -s[1])
        y = 0.0
        for sym, h, col in up:
            _glyph(ax, x, y, sym, h, col, flip=False)
            y += h
        y = 0.0
        for sym, h, col in dn:
            y += h
            _glyph(ax, x, y, sym, -h, col, flip=True)
    n = len(per_col)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(min(-0.2, -ymax * 0.3), ymax * 1.05 or 1)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_ylabel(title, rotation=0, ha='right', va='center', fontsize=9)
    ax.set_yticks([])
    ax.set_xticks(range(0, n, 2))
    ax.tick_params(labelsize=6)


def draw_ic_bars(ax, ic_per_col, title):
    ax.bar(range(len(ic_per_col)), ic_per_col, color='#4c72b0', width=0.8)
    ax.set_ylabel(title, rotation=0, ha='right', va='center', fontsize=9)
    ax.set_xlim(-0.5, len(ic_per_col) - 0.5)
    ax.set_ylim(0, (max(ic_per_col) if ic_per_col else 1) * 1.1 or 1)
    ax.tick_params(labelsize=6)


def draw_lgg_track(ax, labels, title):
    for x, lbl in enumerate(labels):
        specific = '(' in lbl and not lbl.split('(', 1)[1][:1].isupper()
        ax.text(x, 0.02, lbl, rotation=90, ha='center', va='bottom', fontsize=4.4,
                family='monospace', color=('#111' if specific else '#aaa'), clip_on=True)
    ax.set_ylabel(title, rotation=0, ha='right', va='center', fontsize=9)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(labelsize=6)


# --------------------------------------------------------------------------- #
# high-level: build columns + render
# --------------------------------------------------------------------------- #
def _columns(aligned_rows: list[list[int]]) -> list[list[int]]:
    L = max((len(r) for r in aligned_rows), default=0)
    return [[r[p] if p < len(r) else 0 for r in aligned_rows] for p in range(L)]


def relational_logo(
    aligned_rows: list[list[int]],
    atom_store: dict[int, tuple],
    out_path: str,
    *,
    title: str = '',
    levels: list[dict] | None = None,
    color_fn: Callable[[str], str] | None = None,
):
    """Render a relational sequence logo for a multiple alignment.

    Args:
        aligned_rows: aligned int-id rows (gap = 0) from ``progressive_msa``.
        atom_store: id -> structured tuple atom (predicate, *args).
        out_path: output image path (matplotlib figure).
        title: figure title.
        levels: abstraction levels to render as stacked logos. Each is a dict
            ``{'name': str, 'symbol': callable(tuple)->str}``. Default is a single
            "ground" level (the full atom term). Add e.g. a predicate-only level
            for a coarser abstraction.
        color_fn: symbol-string -> colour; default assigns a stable palette.

    Returns the per-column ground information content (list of bits). Also writes
    ``<out_path>`` (figure) and ``<out_path>.ic.csv`` (numeric logo data).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if levels is None:
        levels = [{'name': 'ground', 'symbol': term_str}]
    cmap: dict[str, str] = {}

    def color(sym: str) -> str:
        if color_fn is not None:
            return color_fn(sym)
        if sym not in cmap:
            cmap[sym] = _PALETTE[len(cmap) % len(_PALETTE)]
        return cmap[sym]

    cols = _columns(aligned_rows)
    L = len(cols)

    def atom(i):
        return None if i == 0 else atom_store.get(i)

    # per-level stacked logos
    level_tracks = []
    for lv in levels:
        symfn = lv['symbol']
        syms_per_atom = {i: symfn(a) for i, a in atom_store.items()}
        alpha = len(set(syms_per_atom.values())) or 1
        per_col = []
        for col in cols:
            symbols = [syms_per_atom[i] if i in syms_per_atom else '-' for i in col]
            per_col.append([(s, ic, color(s)) for s, ic in column_ic(symbols, alpha)])
        level_tracks.append((lv['name'], per_col))

    # ground IC bars + per-column lgg
    ground = level_tracks[0][1]
    ground_ic = [sum(h for _, h, _ in ground[p] if h > 0) for p in range(L)]
    lgg_labels = [term_str(lgg_atoms([atom(i) for i in col if i])) for col in cols]

    n_tracks = len(level_tracks) + 2  # + IC bars + lgg
    fig, axes = plt.subplots(n_tracks, 1, figsize=(max(8, L * 0.34), 2 + 1.6 * n_tracks),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, (name, per_col) in zip(axes, level_tracks):
        ymax = max((sum(h for _, h, _ in st if h > 0) for st in per_col), default=1.0) or 1.0
        draw_stack_logo(ax, per_col, f'{name}\n(bits)', ymax)
    draw_ic_bars(axes[len(level_tracks)], ground_ic, 'ground IC\n(bits)')
    draw_lgg_track(axes[-1], lgg_labels, 'lgg\n(per column)')
    axes[-1].set_xlabel('aligned position')
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0.04, 0, 1, 0.98))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)

    csv = f'{out_path}.ic.csv'
    with open(csv, 'w') as fh:
        fh.write('pos,ground_IC_bits,lgg\n')
        for p in range(L):
            fh.write(f'{p},{ground_ic[p]:.4f},"{lgg_labels[p]}"\n')
    return ground_ic
