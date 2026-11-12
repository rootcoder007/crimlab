# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Disparity tests (R/disp.R)."""

from __future__ import annotations

from typing import Any

from . import _clab
from ._out import Table, Vector, tab, vec
from ._util import Vec, grp, items, named, order, rows, twice


def _split2(day: Any) -> tuple[Any, Any]:
    """vodtst(day=<table>) with night left out: one table whose two columns
    are named day and night, in either order (never taken by position): a
    mapping of the two columns (a dict, a pandas DataFrame), each a
    sequence or a mapping from group to stops; or a sequence of rows that
    name them."""
    p = items(day)
    if p is not None and sorted(str(k) for k, _ in p) == ["day", "night"]:
        d = dict(p)
        return d["day"], d["night"]
    r, g, cols = rows(day) if p is None else (None, None, None)
    if r is None or cols is None or sorted(map(str, cols)) != ["day", "night"]:
        raise ValueError("vodtst: `night` is missing: give `night`, or `day` as one table (a "
                         "matrix or a data frame) with two columns named day and night.")
    out = []
    for c in ("day", "night"):
        v = [e[cols.index(c)] for e in r]
        out.append(v if g is None else list(zip(g, v)))
    return out[0], out[1]


def vodtst(day: Vec, night: Vec | None = None, ref: Any = 1,
           shift: Vec | None = None) -> Table:
    """Veil-of-darkness test.

    Compares each group's share of stops made in daylight, when a driver's
    race can be seen, with its share in darkness, when it cannot. The ratio
    of the two odds against a reference group is the factor by which
    visibility changes that group's chance of being stopped; 1 means no
    change. It identifies the visibility effect only when the mix of
    drivers on the road is the same in daylight and darkness: ``shift``
    states how much each group's exposure differs between the two, and
    ``adj`` is the ratio with that difference taken out. See
    ``proof("vodtst")``.

    Parameters
    ----------
    day : stops by group in daylight, non-negative, optionally named.
        With ``night`` left out, one table of groups by stops whose two
        columns are named ``day`` and ``night``: a mapping of the two
        columns, each a sequence or a mapping from group to stops.
    night : stops by group in darkness, positive, one per group; matched to
        ``day`` by name when both are named.
    ref : the reference group, by name (matched exactly) or by position
        from 1 in ``day`` as given.
    shift : None, or one positive number per group: daylight exposure over
        darkness exposure.

    Returns
    -------
    A data frame (dict of columns), one row per group in the order of
    ``day``: ``group``, ``ref`` (True for the reference group), ``shday``,
    ``shnight``, ``ratio``, ``sign`` (the sign of ``shday - shnight``)
    and, with ``shift``, ``adj``.

    Examples
    --------
    >>> r = vodtst(day={"white": 600, "black": 400},
    ...            night={"white": 700, "black": 300})
    >>> r["group"], [round(x, 4) for x in r["ratio"]]
    (['white', 'black'], [1.0, 1.5556])
    """
    if night is None:
        day, night = _split2(day)
    (d, n, s), g = named("vodtst", [("day", day), ("night", night), ("shift", shift)])
    return tab(_clab.vodtst(d, n, ref, s, g))


def orrisk(*args: float, **kwargs: float) -> Vector:
    """Risks behind an odds ratio.

    An odds ratio is not a relative risk: it lies further from 1, and by
    more the more common the outcome. Given the odds ratio, the overall
    rate of the outcome and the share exposed, this returns the two risks
    and the relative risk they imply. See ``proof("orrisk")``.

    The arguments are ``or``, ``rate`` and ``share``, as in R. ``or`` is a
    Python keyword, so it is given by position or as ``**{"or": 3}``.

    Parameters
    ----------
    or : the odds ratio, one positive number.
    rate : the overall rate of the outcome, strictly between 0 and 1.
    share : the exposed share of the population, strictly between 0 and 1.

    Returns
    -------
    ``exposed`` and ``unexposed`` (the two risks), ``rr`` (their ratio) and
    ``or``. A risk closer to 1 than a double can hold is returned as 1.

    Examples
    --------
    >>> {k: round(v, 4) for k, v in orrisk(3, 0.4, 0.5).items()}
    {'exposed': 0.5282, 'unexposed': 0.2718, 'rr': 1.9436, 'or': 3.0}
    """
    names = ("or", "rate", "share")
    extra = [k for k in kwargs if k not in names]
    if len(args) > 3 or extra:
        raise TypeError(f"orrisk: takes or, rate and share, not {extra or len(args)}")
    a = dict(zip(names, args))
    for k, v in kwargs.items():
        if k in a:
            raise TypeError(f"orrisk: `{k}` given twice")
        a[k] = v
    miss = [k for k in names if k not in a]
    if miss:
        raise TypeError(f"orrisk: missing {', '.join(miss)}")
    return vec(_clab.orrisk(a["or"], a["rate"], a["share"]))


def outtst(q: Vec, t: float, w: Vec | None = None) -> Vector:
    """Outcome test.

    The stop rate and the hit rate of a group searched whenever the chance
    of guilt exceeds a threshold. Two groups held to the same threshold can
    show different hit rates (infra-marginality). See ``proof("outtst")``.

    Parameters
    ----------
    q : the chances of guilt in the group, from 0 to 1.
    t : the search threshold, one number.
    w : None for equal weights, or one non-negative weight per element of
        ``q`` with a positive total; matched to ``q`` by name when both are
        named with unique names.

    Returns
    -------
    ``stop`` (the weighted share searched) and ``hit`` (the weighted mean
    chance of guilt among them; NaN when nobody is searched).

    Examples
    --------
    >>> outtst([0.1, 0.5, 0.9], 0.3)
         stop       hit
    0.6666667 0.7000000
    """
    (qv, wv), _ = named("outtst", [("q", q), ("w", w)], strict=False)
    return vec(_clab.outtst(qv, t, wv))


def dspbm(pop: Vec, contact: Vec, force: Vec, ref: Any = 1) -> Table:
    """Benchmark decomposition of a disparity.

    Splits each group's disparity in force per resident, against a
    reference group, into a disparity in contact per resident times a
    disparity in force per contact; the split is a product, never a sum.
    See ``proof("dspbm")``.

    Parameters
    ----------
    pop : residents by group, positive, optionally named.
    contact : police contacts by group, positive, one per group.
    force : uses of force by group, positive, one per group.
    ref : the reference group, by name or by position from 1.

    Returns
    -------
    A data frame in the order of ``pop``: ``group``, ``ref`` (True for the
    reference group), ``cd`` (contact per resident), ``fd`` (force
    per contact) and ``rd`` (force per resident), each against the
    reference group, with ``rd = cd * fd``.

    Examples
    --------
    >>> r = dspbm(pop={"a": 9000, "b": 1000}, contact={"a": 900, "b": 300},
    ...           force={"a": 9, "b": 6})
    >>> [round(x, 4) for x in r["rd"]]
    [1.0, 6.0]
    """
    (p, c, f), g = named("dspbm", [("pop", pop), ("contact", contact), ("force", force)])
    return tab(_clab.dspbm(p, c, f, ref, g))


def dyadrt(count: Any, pop: Vec) -> Table:
    """Dyadic offending rates.

    Offences by one group against another, as a rate per member of the
    offending group and as a rate per offender-victim pair exposed. The
    rate per offender rises with the size of the victim group even when
    offenders choose victims at random; the rate per pair does not. See
    ``proof("dyadrt")``.

    Parameters
    ----------
    count : a square table of offence counts, offending group in rows and
        victim group in columns: a sequence of rows, a 2-D array, or a
        mapping of columns ``{victim group: {offending group: count}}``. When it
        names both its rows and its columns, the columns are matched to the
        rows by name, and both to ``pop`` when ``pop`` is named.
    pop : group sizes, positive, optionally named; unnamed, in the order of
        the rows.

    Returns
    -------
    A data frame with one row per dyad: ``off``, ``vic``, ``count``,
    ``rtoff``, ``rtpair``, ``nulloff`` (the rate per offender under random
    mixing) and ``ratio`` (``rtpair`` over the pooled pair rate).

    Examples
    --------
    >>> r = dyadrt([[640, 160], [160, 40]], {"a": 800, "b": 200})
    >>> r["off"], r["vic"], [round(x, 4) for x in r["ratio"]]
    (['a', 'b', 'a', 'b'], ['a', 'a', 'b', 'b'], [1.0, 1.0, 1.0, 1.0])
    """
    p, pn = grp(pop)
    tab, rn, cn = rows(count)
    rd, cd = "the row names of `count`", "the column names of `count`"
    for what, nm in (("`pop`", pn), (rd, rn), (cd, cn)):
        k = twice(nm or [])
        if k is not None:
            raise ValueError(f"dyadrt: {what} names a group twice: {k}.")
    g = len(p) if hasattr(p, "__len__") else 1
    if tab is not None and len(tab) == g and all(len(r) == g for r in tab) \
            and rn is not None and cn is not None:
        # the columns in the order of the rows, both in the order of pop
        i = list(range(g))
        ro = order("dyadrt", "`pop`", pn, rd, rn, i) if pn is not None else i
        co = order("dyadrt", "`pop`", pn, cd, cn, i) if pn is not None else \
            order("dyadrt", rd, rn, cd, cn, i)
        tab = [[tab[r][c] for c in co] for r in ro]
    names = pn if pn is not None else rn if rn is not None else cn
    nr = len(tab) if tab is not None else -1
    nc = len(tab[0]) if tab else (0 if tab is not None else -1)
    if tab is not None and any(len(r) != nc for r in tab):
        tab = None
    flat = None if tab is None else [tab[i][j] for j in range(nc) for i in range(nr)]
    return Table(_clab.dyadrt(flat, nr, nc, p, names))
