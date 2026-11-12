# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Forensic statistics (R/frns.R, R/frns2.R)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from . import _clab
from ._util import Frame, Named, Vec, align, flat, isseq, named, nrows, rows, table, vals

#: The prior log10 odds eceplt() draws at by default: R's
#: seq(-2.5, 2.5, by = 0.05).
LO = [min(-2.5 + i * 0.05, 2.5) for i in range(101)]


def gtmp(p: Vec, q: Vec | None = None, hom: Any = False, theta: float = 0) -> dict[str, Any]:
    """Genotype match probability.

    The chance that an unrelated person from the same population shares a
    DNA profile, locus by locus and overall. With ``theta = 0`` this is the
    product rule (``p^2`` for a homozygote, ``2pq`` for a heterozygote); a
    positive ``theta`` applies the Balding-Nichols correction for shared
    ancestry. See ``proof("gtmp")``.

    Parameters
    ----------
    p : frequency of the first allele at each locus, above 0 and at most 1,
        optionally named by locus.
    q : frequency of the second allele at each locus; not read where
        ``hom`` is true (None there is fine).
    hom : True where the genotype is homozygous: one bool, or one per locus.
    theta : the coancestry coefficient, from 0 up to, not including, 1.

    Returns
    -------
    ``locus``, a data frame with ``locus``, ``hom`` and ``mp``;
    ``log10mp``, the base-10 logarithm of the profile match probability;
    and ``log10lr``, the same with the sign changed.

    Examples
    --------
    >>> r = gtmp(p={"D3": 0.1, "vWA": 0.2}, q=[0.3, None], hom=[False, True],
    ...          theta=0.03)
    >>> r["locus"]["locus"], round(r["log10lr"], 4)
    (['D3', 'vWA'], 2.3068)
    """
    (pv, qv, hv), g = named("gtmp", [("p", p), ("q", q), ("hom", hom)])
    return _clab.gtmp(pv, qv, hv, vals(theta), g)


def cllr(same: Vec, diff: Vec) -> float:
    """Log-likelihood-ratio cost.

    Scores a set of likelihood ratios against the truth: 1 for a system
    that always answers 1, above 0 for any set of ratios, tending to 0 for
    a system that separates the two kinds perfectly. See ``proof("cllr")``
    and :func:`cllrmin`.

    Parameters
    ----------
    same : likelihood ratios from same-source comparisons, positive.
    diff : likelihood ratios from different-source comparisons, positive.

    Returns
    -------
    The cost in bits.

    Examples
    --------
    >>> round(cllr([50, 200, 8], [0.01, 0.2, 1.5]), 6)
    0.300835
    """
    return _clab.cllr(vals(same), vals(diff))


def lrpost(lr: Vec, prior: Vec) -> list[float]:
    """Posterior probability from a likelihood ratio.

    A large ratio with a small prior can still leave the probability of the
    first hypothesis low; reading the ratio as that probability is the
    prosecutor's fallacy. See ``proof("lrpost")``.

    Parameters
    ----------
    lr : likelihood ratios, positive.
    prior : the prior probability of the first hypothesis, strictly between
        0 and 1: one value, or one per ratio, matched to ``lr`` by name when
        both are named with unique names.

    Returns
    -------
    The posterior probabilities, one per ratio.

    Examples
    --------
    >>> [round(x, 6) for x in lrpost(1e6, 1 / 5e6)]
    [0.166667]
    """
    (lv, pv), _ = named("lrpost", [("lr", lr), ("prior", prior)], strict=False)
    return _clab.lrpost(lv, pv)


def _names(x: Any) -> Any:
    """Allele names: one text is one name; any other collection a list."""
    if isinstance(x, str):
        return [x]
    x = flat(x)
    if isseq(x) or isinstance(x, (set, frozenset)):
        return list(x)
    return x


def mixlr(p: Mapping[Any, float] | Any, np: str | Sequence[Any], xp: int,
          nd: str | Sequence[Any], xd: int) -> Named:
    """Mixture likelihood ratio.

    The weight of a DNA mixture at one locus for one account of who
    contributed against another. Each account names the alleles of the
    mixture its known contributors do not explain and the number of
    unknown contributors who must explain them. Peak heights, drop-out and
    coancestry are not modelled. See ``proof("mixlr")``.

    Parameters
    ----------
    p : the population frequencies of the alleles in the mixture, named
        (a mapping or a Series), above 0 with a sum of at most 1.
    np, nd : the alleles left unexplained under the first and the second
        account: names of ``p`` (one name, or a list, possibly empty), at
        most 20.
    xp, xd : the number of unknown contributors under each account.

    Returns
    -------
    ``num``, ``den`` and ``lr``.

    Examples
    --------
    >>> p = {"a": 0.1, "b": 0.2, "c": 0.05, "d": 0.15}
    >>> {k: round(v, 4) for k, v in mixlr(p, [], 0, ["c", "d"], 1).items()}
    {'num': 1.0, 'den': 0.015, 'lr': 66.6667}
    """
    (pv,), g = named("mixlr", [("p", p)])
    return _clab.mixlr(pv, g, _names(np), vals(xp), _names(nd), vals(xd))


def dbsrch(mp: float, n: float, pop: float | None = None) -> Named:
    """Database search.

    What a match found by trawling a database is worth: ``np``, the match
    probability times the number of profiles searched (the second National
    Research Council report); ``any``, the exact chance of at least one
    chance match, never above ``np``; and ``post``, the matching person's
    probability of being the source among ``pop`` possible sources with
    equal prior (Balding and Donnelly). See ``proof("dbsrch")``.

    Parameters
    ----------
    mp : the profile match probability, above 0 and at most 1.
    n : the number of profiles searched, the match included.
    pop : None, or the number of possible sources, at least ``n``.

    Returns
    -------
    ``np``, ``any`` and, with ``pop``, ``post``.

    Examples
    --------
    >>> {k: round(v, 6) for k, v in dbsrch(1e-9, 5e6, pop=3e7).items()}
    {'np': 0.005, 'any': 0.004988, 'post': 0.97561}
    """
    return _clab.dbsrch(vals(mp), vals(n), vals(pop))


def cllrmin(same: Vec, diff: Vec) -> Named:
    """Log-likelihood-ratio cost after recalibration.

    Splits the cost of a set of likelihood ratios into ``min``, their cost
    after recalibration by the pool-adjacent-violators algorithm, and
    ``cal``, the rest (``cllr - min``). See ``proof("cllrmin")`` and
    ``proof("lrcal")``.

    Parameters
    ----------
    same, diff : as for :func:`cllr`.

    Returns
    -------
    ``cllr``, ``min`` and ``cal``, with ``cllr = min + cal``.

    Examples
    --------
    >>> {k: round(v, 6) for k, v in cllrmin([50, 200, 8], [0.01, 0.2, 1.5]).items()}
    {'cllr': 0.300835, 'min': 0.0, 'cal': 0.300835}
    """
    return _clab.cllrmin(vals(same), vals(diff))


def lrcal(same: Vec, diff: Vec) -> dict[str, list[float]]:
    """Calibrated likelihood ratios.

    Recalibrates a set of likelihood ratios against the truth with the
    pool-adjacent-violators algorithm, keeping their order. The result can
    hold 0 and infinity: the data never contradicted those ratios. See
    ``proof("lrcal")``.

    Parameters
    ----------
    same, diff : as for :func:`cllr`.

    Returns
    -------
    ``same`` and ``diff``, the calibrated ratios in the order given.

    Examples
    --------
    >>> lrcal([3, 0.7, 20, 5], [0.2, 1, 4, 0.05])
    {'same': [1.0, 1.0, inf, inf], 'diff': [0.0, 1.0, 1.0, 0.0]}
    """
    return _clab.lrcal(vals(same), vals(diff))


def _axes(draw: bool, ax: Any) -> Any:
    """The axes to draw on: ax, or the current ones; None when there is
    nothing to draw or matplotlib is not installed."""
    if not draw:
        return None
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt.gca()


def tipplt(same: Vec, diff: Vec, *, draw: bool = True, ax: Any = None,
           **kwargs: Any) -> dict[str, Frame]:
    """Tippett plot.

    For each value on the axis, the share of same-source and of
    different-source likelihood ratios at or above it. ``proof()`` records
    no theorem for this plot.

    Parameters
    ----------
    same, diff : as for :func:`cllr`.
    draw : draw the plot (when matplotlib is installed).
    ax : the matplotlib axes to draw on; the current axes by default.
    **kwargs : passed to ``ax.set()`` (``title=...``, ``xlim=...``).

    Returns
    -------
    ``same`` and ``diff``, each a data frame with ``log10lr`` and ``share``
    (the share of ratios at or above it).

    Examples
    --------
    >>> r = tipplt([50, 200, 8, 0.5], [0.01, 0.2, 1.5, 0.001], draw=False)
    >>> r["same"]["share"]
    [1.0, 0.75, 0.5, 0.25]
    """
    r = _clab.tipplt(vals(same), vals(diff))
    a = _axes(draw, ax)
    if a is not None:
        s, d = r["same"], r["diff"]
        a.step(s["log10lr"], s["share"], where="post", linestyle="-", label="same source")
        a.step(d["log10lr"], d["share"], where="post", linestyle="--",
               label="different source")
        a.axvline(0, linestyle=":")
        a.set(xlabel="log10 likelihood ratio", ylabel="share of ratios at or above",
              ylim=(0, 1), **kwargs)
        a.legend(loc="lower left", frameon=False)
    return r


def eceplt(same: Vec, diff: Vec, lo: Vec = LO, *, draw: bool = True, ax: Any = None,
           **kwargs: Any) -> Frame:
    """Empirical cross-entropy plot.

    The information a set of likelihood ratios still leaves missing, across
    prior odds: for the ratios as reported, after calibration, and for a
    system that always reports 1. At even prior odds the first two curves
    pass through the values of :func:`cllrmin`. ``proof()`` records no
    theorem for this plot.

    Parameters
    ----------
    same, diff : as for :func:`cllr`.
    lo : the prior log10 odds to draw at; by default -2.5 to 2.5 by 0.05.
    draw, ax, **kwargs : as for :func:`tipplt`.

    Returns
    -------
    A data frame: ``lo``, ``ece``, ``min`` and ``ref``.

    Examples
    --------
    >>> r = eceplt([50, 200, 8, 0.5], [0.01, 0.2, 1.5, 0.001], lo=[0], draw=False)
    >>> [round(r[k][0], 6) for k in ("ece", "min", "ref")]
    [0.423926, 0.25, 1.0]
    """
    r = _clab.eceplt(vals(same), vals(diff), vals(lo))
    a = _axes(draw, ax)
    if a is not None:
        a.plot(r["lo"], r["ece"], linestyle="-", label="reported")
        a.plot(r["lo"], r["min"], linestyle="--", label="calibrated")
        a.plot(r["lo"], r["ref"], linestyle=":", label="always 1")
        a.set(xlabel="prior log10 odds", ylabel="empirical cross-entropy (bits)",
              ylim=(0, max(r["ece"] + r["ref"])), **kwargs)
        a.legend(loc="upper left", frameon=False)
    return r


def _freq(f: Any) -> tuple[list[Any] | None, list[Any] | None]:
    """freq as rows of two, and its row names; rows None when it is not a
    table of pairs."""
    r, fl = rows(f)
    if r is None or any(len(e) != 2 for e in r):
        return None, fl
    return r, fl


def _codes(*tabs: list[str]) -> list[list[int]]:
    """The alleles of the tables as codes, one per distinct label."""
    u: dict[str, int] = {}
    return [[u.setdefault(a, len(u) + 1) for a in t] for t in tabs]


def _pairs(v: list[Any]) -> list[list[Any]]:
    return [v[i:i + 2] for i in range(0, len(v), 2)]


def _loci(fn: str, tabs: list[tuple[str, list[str], list[Any] | None]], freq: Any,
          gn: str) -> tuple[list[list[str]], list[Any] | None, list[Any] | None]:
    """The tables (allele labels, two per row) and freq matched by locus
    under R's rule (align()): the tables, freq row by row (None when it is
    no table of pairs) and the loci, named by the first table (freq
    included) that names them."""
    fr, fl = _freq(freq)
    args = [(nm, _pairs(t), loc) for nm, t, loc in tabs] + [("freq", fr, fl)]
    vs, loci = align(fn, args, True)
    out = [[a for e in v for a in e] for v in vs[:-1]]
    f = None if vs[-1] is None else [a for e in vs[-1] for a in e]
    return out, f, loci


def kinidx(g1: Any, g2: Any, freq: Any, rel: Any = "parent") -> dict[str, Any]:
    """Kinship index.

    The weight of two DNA profiles for a stated relationship against no
    relationship, locus by locus and overall. Loci are taken as
    independent; there is no allowance for mutation or coancestry. See
    ``proof("kinidx")``.

    Parameters
    ----------
    g1, g2 : the two genotypes: tables of allele labels with two columns,
        one row per locus (a sequence of pairs, a 2-D array, a DataFrame,
        or a mapping from locus to pair), or two alleles for one locus.
        When both name their loci, ``g2`` is matched to ``g1`` by name.
    freq : the population frequencies of the alleles of ``g1``, same shape.
    rel : ``"parent"``, ``"sibling"``, ``"half"`` (half siblings,
        grandparent or uncle) or ``"cousin"``, or the three probabilities of
        sharing 0, 1 and 2 alleles by descent.

    Returns
    -------
    ``locus``, a data frame with ``locus`` and ``idx``; ``log10idx``, the
    base-10 logarithm of the combined index; and ``post``, the probability
    of the relationship at even prior odds.

    Examples
    --------
    >>> r = kinidx([12, 14], [12, 14], [0.1, 0.2], rel="sibling")
    >>> [round(x, 4) for x in r["locus"]["idx"]]
    [8.375]
    """
    a, l1 = table("kinidx", "g1", g1)
    b, l2 = table("kinidx", "g2", g2)
    nrows("kinidx", "g2", b, len(a) // 2)
    (a, b), f, loci = _loci("kinidx", [("g1", a, l1), ("g2", b, l2)], freq, "g1")
    c = _codes(a, b)
    return _clab.kinidx(c[0], c[1], f, flat(rel), loci)


def patidx(child: Any, father: Any, freq: Any, mother: Any = None) -> dict[str, Any]:
    """Paternity index.

    The weight of the DNA evidence that an alleged father is the father of
    a child, against an unrelated man, locus by locus and overall. With the
    mother's genotype the index uses the allele the child must have
    received from the father; without it the index is the parent-child
    :func:`kinidx`. One inconsistent locus gives an index of 0 (no
    allowance for mutation). See ``proof("patidx")``.

    Parameters
    ----------
    child, father : genotypes, as for :func:`kinidx`.
    freq : the population frequencies of the child's alleles, same shape.
    mother : None, or the mother's genotype in the same form.

    Returns
    -------
    ``locus``, a data frame with ``locus`` and ``idx`` (NaN where the
    mother's genotype rules her out as mother); ``log10idx``; and ``post``,
    the probability of paternity at even prior odds.

    Examples
    --------
    >>> r = patidx([12, 14], [14, 15], [0.1, 0.2], mother=[12, 12])
    >>> [round(x, 4) for x in r["locus"]["idx"]]
    [2.5]
    """
    if mother is None:
        # kinidx(child, father, freq, "parent"): its messages name g1, g2
        a, l1 = table("patidx", "g1", child)
        b, l2 = table("patidx", "g2", father)
        nrows("patidx", "g2", b, len(a) // 2)
        (a, b), f, loci = _loci("patidx", [("g1", a, l1), ("g2", b, l2)], freq, "g1")
        c = _codes(a, b)
        return _clab.patidx(c[0], c[1], None, f, loci)
    a, l1 = table("patidx", "child", child)
    b, l2 = table("patidx", "father", father)
    nrows("patidx", "father", b, len(a) // 2)
    m, l3 = table("patidx", "mother", mother)
    nrows("patidx", "mother", m, len(a) // 2)
    (a, b, m), f, loci = _loci("patidx", [("child", a, l1), ("father", b, l2),
                                          ("mother", m, l3)], freq, "child")
    c = _codes(a, b, m)
    return _clab.patidx(c[0], c[1], c[2], f, loci)
