"""Frequency and probability distributions for relational sequence logos.

Simplified, modern-Python reimplementation of the NLTK-derived legacy
``Probability.py``.  Only the distributions needed for logo construction
are included: :class:`FreqDist`, :class:`MLEProbDist` (maximum-likelihood),
and :class:`LidstoneProbDist` (smoothed).
"""

from __future__ import annotations

import builtins as _builtins
import math
from collections.abc import Hashable, Iterator
from typing import Any

_builtin_max = _builtins.max


class FreqDist:
    """A frequency distribution over hashable samples.

    Counts how many times each outcome has been observed.

    >>> fd = FreqDist()
    >>> fd.inc("a"); fd.inc("a"); fd.inc("b")
    >>> fd.count("a")
    2
    >>> fd.freq("a")  # doctest: +ELLIPSIS
    0.666...
    """

    __slots__ = ('_counts', '_total', '_max_cache')

    def __init__(self) -> None:
        self._counts: dict[Any, int] = {}
        self._total: int = 0
        self._max_cache: Any | None = None

    # -- Mutation -----------------------------------------------------------

    def inc(self, sample: Hashable, count: int = 1) -> None:
        """Increment the count for *sample* by *count*."""
        if count == 0:
            return
        self._counts[sample] = self._counts.get(sample, 0) + count
        self._total += count
        self._max_cache = None

    # -- Queries ------------------------------------------------------------

    @property
    def total(self) -> int:
        """Total number of recorded outcomes (``N``)."""
        return self._total

    @property
    def num_bins(self) -> int:
        """Number of distinct samples with count > 0 (``B``)."""
        return len(self._counts)

    def count(self, sample: Hashable) -> int:
        """Return the count for *sample* (0 if unseen)."""
        return self._counts.get(sample, 0)

    def freq(self, sample: Hashable) -> float:
        """Return the relative frequency ``count(sample) / N``."""
        if self._total == 0:
            return 0.0
        return self._counts.get(sample, 0) / self._total

    def samples(self) -> list[Any]:
        """Return all samples with count > 0."""
        return list(self._counts.keys())

    def max(self) -> Any | None:
        """Return the sample with the highest count (arbitrary tie-break)."""
        if self._max_cache is None:
            if not self._counts:
                return None
            self._max_cache = _builtin_max(self._counts, key=self._counts.__getitem__)
        return self._max_cache

    def sorted_samples(self) -> list[Any]:
        """Return samples sorted by descending count."""
        return sorted(self._counts, key=self._counts.__getitem__, reverse=True)

    # -- Container protocol -------------------------------------------------

    def __contains__(self, sample: object) -> bool:
        return sample in self._counts

    def __len__(self) -> int:
        return self.num_bins

    def __iter__(self) -> Iterator[Any]:
        return iter(self._counts)

    # -- Representation -----------------------------------------------------

    def __repr__(self) -> str:
        return f'<FreqDist with {self._total} outcomes, {self.num_bins} bins>'

    def __str__(self) -> str:
        items = ', '.join(
            f'{s!r}: {c}' for s, c in sorted(self._counts.items(), key=lambda kv: -kv[1])
        )
        return f'<FreqDist: {items}>'


class MLEProbDist:
    """Maximum-likelihood probability distribution from a :class:`FreqDist`.

    ``P(sample) = count(sample) / N``
    """

    __slots__ = ('_freqdist',)

    def __init__(self, freqdist: FreqDist) -> None:
        if freqdist.total == 0:
            raise ValueError('Cannot build MLE distribution from empty FreqDist.')
        self._freqdist = freqdist

    @property
    def freqdist(self) -> FreqDist:
        return self._freqdist

    def prob(self, sample: Hashable) -> float:
        return self._freqdist.freq(sample)

    def logprob(self, sample: Hashable) -> float:
        p = self.prob(sample)
        return math.log(p) if p > 0 else float('-inf')

    def max(self) -> Any | None:
        return self._freqdist.max()

    def samples(self) -> list[Any]:
        return self._freqdist.samples()

    def __repr__(self) -> str:
        return f'<MLEProbDist based on {self._freqdist.total} outcomes>'


class LidstoneProbDist:
    """Lidstone-smoothed probability distribution.

    ``P(sample) = (count(sample) + gamma) / (N + B * gamma)``

    With ``gamma = 1`` this is Laplace smoothing; ``gamma = 0.5`` gives the
    Expected Likelihood Estimate (ELE).
    """

    __slots__ = ('_freqdist', '_gamma', '_bins', '_N')

    def __init__(
        self,
        freqdist: FreqDist,
        gamma: float = 1.0,
        bins: int | None = None,
    ) -> None:
        if bins is not None and bins < freqdist.num_bins:
            raise ValueError(f'bins ({bins}) must be >= FreqDist.num_bins ({freqdist.num_bins})')
        if bins is None:
            bins = freqdist.num_bins
        if bins == 0:
            raise ValueError('Lidstone distribution must have at least one bin.')

        self._freqdist = freqdist
        self._gamma = float(gamma)
        self._bins = bins
        self._N = freqdist.total

    @property
    def freqdist(self) -> FreqDist:
        return self._freqdist

    def prob(self, sample: Hashable) -> float:
        c = self._freqdist.count(sample)
        return (c + self._gamma) / (self._N + self._bins * self._gamma)

    def logprob(self, sample: Hashable) -> float:
        p = self.prob(sample)
        return math.log(p) if p > 0 else float('-inf')

    def max(self) -> Any | None:
        return self._freqdist.max()

    def samples(self) -> list[Any]:
        return self._freqdist.samples()

    def __repr__(self) -> str:
        return f'<LidstoneProbDist gamma={self._gamma} based on {self._N} outcomes>'
