# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Corrections and sentencing (R/corr.R)."""

from __future__ import annotations

from typing import Any

from . import _clab
from ._out import Table, Vector, tab, vec
from ._util import Vec, grp, label1, labels, named


def incap(lam: Vec, q: Vec, s: Vec, share: Vec | None = None) -> Table:
    """Incapacitation: crime rate, share prevented and the marginal year.

    An offender commits crimes at rate ``lam`` a year while free, and each
    crime leads to a custodial sentence of ``s`` years with probability
    ``q``. Over a career the crime rate is ``lam / (1 + lam q s)``; the
    share prevented, ``lam q s / (1 + lam q s)``, never reaches 1, and one
    more year of sentence prevents less for each year added. See
    ``proof("incap")`` and :func:`incapr`.

    Parameters
    ----------
    lam : offending rates while free (crimes a year), positive, one per
        group, optionally named.
    q : the chance a crime leads to a custodial sentence, above 0 and at
        most 1: one value, or one per group.
    s : sentence lengths in years, non-negative: one value, or one per
        group.
    share : None, or one non-negative population share per group, summing
        to 1.

    Returns
    -------
    A data frame: ``group``, ``lam``, ``q``, ``s``, ``rate`` (the long-run
    crime rate), ``prev`` (the share prevented) and ``marg`` (the crimes a
    year one more year of sentence prevents). With ``share``, a last row
    ``all`` holds the weighted free rate in ``lam``, the weighted rate in
    ``rate`` and the share of all offending prevented in ``prev``; its
    ``q``, ``s`` and ``marg`` are None (R's NA).

    Examples
    --------
    >>> r = incap(lam={"low": 2, "high": 10}, q=0.1, s=1, share=[0.8, 0.2])
    >>> r["group"], [round(x, 4) for x in r["prev"]]
    (['low', 'high', 'all'], [0.1667, 0.5, 0.3519])
    """
    (lv, qv, sv, hv), g = named("incap", [("lam", lam), ("q", q), ("s", s), ("share", share)])
    return tab(_clab.incap(lv, qv, sv, hv, g))


def incapr(lam: Vec, t0: int, s: int, r: float = 0) -> Vector:
    """Incapacitation under desistance and replacement.

    Crimes prevented by a sentence of ``s`` periods served from career age
    ``t0``, when the offending rate falls with age and a share ``r`` of the
    crimes prevented is committed by someone else instead. The constant-rate
    estimate overstates what a sentence buys under both. See
    ``proof("incapr")``.

    Parameters
    ----------
    lam : offending rates by career age, ``lam[0]`` at age 0, non-negative
        and never increasing.
    t0 : the career age at which the sentence starts, a whole number from
        0, with ``t0 + s`` at most ``len(lam) - 1``.
    s : the sentence in periods, a whole number from 1.
    r : the replacement share, from 0 to 1.

    Returns
    -------
    ``prev`` (crimes prevented), ``net`` (net of replacement), ``upper``
    (``s * lam`` at age ``t0``), ``lower`` (``s * lam`` at age ``t0 + s``)
    and ``later`` (crimes prevented by the same sentence from age
    ``t0 + 1``).

    Examples
    --------
    >>> incapr([12, 10, 8, 6, 5, 4, 3, 2, 2, 1], t0=2, s=3, r=0.25)
     prev   net upper lower later
    19.00 14.25 24.00 12.00 15.00
    """
    return vec(_clab.incapr(lam, t0, s, r))


def sellab(y: Vec, rel: Any, rule: Any, w: Vec | None = None) -> Vector:
    """Evaluating a release rule from the released alone.

    Failures are seen only for the defendants a judge released. A rule
    that releases people the judge detained has a failure rate known only
    up to an interval whose width is the share of the rule's releases the
    judge detained. See ``proof("sellab")``.

    Parameters
    ----------
    y : failures, 1 or 0 for each released defendant; values for detained
        defendants are ignored and may be None.
    rel : released by the judge, True or False for each defendant.
    rule : released by the rule, True or False, same length.
    w : None for equal weights, or one non-negative weight per defendant.
        When ``y``, ``rel``, ``rule`` (and ``w``) are all named with unique
        names, they are matched by name.

    Returns
    -------
    ``obs`` (the failure rate among the judge's releases), ``unobs`` (the
    share of the rule's releases the judge detained), ``lower`` and
    ``upper`` (the bounds on the rule's failure rate), ``width`` and
    ``naive`` (the failure rate on the rule's observed releases; NaN when
    none were).

    Examples
    --------
    >>> sellab([1, 0, None, 0], [True, True, False, True],
    ...        [True, False, True, True])
          obs     unobs     lower     upper     width     naive
    0.3333333 0.3333333 0.3333333 0.6666667 0.3333333 0.5000000
    """
    (a, b, c, d), _ = named("sellab", [("y", y), ("rel", rel), ("rule", rule), ("w", w)],
                            strict=False)
    return vec(_clab.sellab(a, b, c, d))


def losbnd(done: Vec, open: Vec) -> Vector:
    """Mean time in custody from completed and current stays.

    The mean length of the completed stays leaves out everyone still
    inside; their elapsed time is a lower bound on their stay, so the mean
    over everyone is at least ``lower``. No upper bound exists without a
    cap on how long the current stays will last. See ``proof("losbnd")``.

    Parameters
    ----------
    done : completed stays, non-negative.
    open : elapsed time of the stays still open, non-negative.

    Returns
    -------
    ``mean`` (of the completed stays), ``age`` (the mean elapsed time of
    the open ones), ``lower`` (the lower bound on the mean over all) and
    ``bias`` (``lower - mean``).

    Examples
    --------
    >>> losbnd([30, 45, 60, 90, 120], [100, 150, 200])
      mean     age   lower    bias
    69.000 150.000  99.375  30.375
    """
    return vec(_clab.losbnd(done, open))


def sntbnd(y: Vec, z: Any, treat: Any = None, w: Vec | None = None) -> Table:
    """Bounds on a sentence's effect without assumptions.

    Each person receives one of two sentences and shows the outcome under
    that one only. With no assumption about who receives which, the effect
    of one sentence against the other lies in an interval of width exactly
    1 that always contains 0. See ``proof("sntbnd")``.

    Parameters
    ----------
    y : outcomes, 0 or 1 for each person.
    z : sentences, exactly two distinct labels, one per person: text,
        numbers, bools or a categorical (read by its labels). Labels are
        written as R's ``as.character()`` writes them (``2.0`` is ``"2"``).
    treat : the label of ``z`` taken as the treatment, matched exactly.
        None is allowed only where the treatment is explicit in ``z``:
        True for a ``z`` of bools, 1 for a ``z`` of 0 and 1; otherwise it
        must be given, since a treatment chosen by sorting would be
        arbitrary. The comparison comes first, the two labels in Unicode
        code-point order otherwise.
    w : None for equal weights, or one non-negative weight per person with
        a positive total.

    Returns
    -------
    A data frame with three rows, the comparison sentence, the treatment and
    ``effect``: ``term``, ``share`` (None for the effect), ``lower``,
    ``upper`` and ``observed``.

    Examples
    --------
    >>> r = sntbnd([1, 0, 1, 1, 0, 0], ["a", "a", "b", "b", "b", "a"], treat="b")
    >>> r["term"], [round(x, 4) for x in r["lower"]]
    (['a', 'b', 'effect'], [0.1667, 0.3333, -0.3333])
    """
    (a, b, c), _ = named("sntbnd", [("y", y), ("z", z), ("w", w)], strict=False)
    lv, kind = labels("sntbnd", "z", b)
    return tab(_clab.sntbnd(a, lv, None if treat is None else label1(treat), c, kind))


def hzdsel(s: float, h: float, l: float, surh: Vec, surl: Vec) -> Vector:
    """Period-two hazards after selective depletion.

    A share ``s`` of a cohort is at high risk, ``h`` a period, and the rest
    at low risk, ``l``. If a programme changes who survives the first
    period, the period-two hazard among the survivors rises with the share
    of high-risk survivors: a period-by-period hazard ratio is not a causal
    contrast. See ``proof("hzdsel")``.

    Parameters
    ----------
    s : the high-risk share at the start, strictly between 0 and 1.
    h : the period hazard of the high-risk type, at most 1.
    l : the period hazard of the low-risk type, from 0, below ``h``.
    surh : period-one survival of the high-risk type in the control and
        treated arms, two numbers above 0 and at most 1, the control arm
        first, or named ``control`` and ``treated`` in either order (with
        ``surl`` named too).
    surl : period-one survival of the low-risk type, the same way; named
        inputs are matched to ``surh`` by name.

    Returns
    -------
    ``share0``, ``share1`` (the high-risk share among the survivors of each
    arm), ``haz0``, ``haz1`` (the period-two hazard in each arm) and
    ``ratio`` (``haz1 / haz0``).

    Examples
    --------
    >>> {k: round(v, 4) for k, v in hzdsel(0.5, 0.5, 0.1, [0.5, 0.8], [0.9, 0.9]).items()}
    {'share0': 0.3571, 'share1': 0.4706, 'haz0': 0.2429, 'haz1': 0.2882, 'ratio': 1.1869}
    """
    (a, b), g = named("hzdsel", [("surh", surh), ("surl", surl)])
    # arms named control and treated (both arguments named) are read by
    # name, in either order; otherwise the first is the control arm
    if g == ["treated", "control"] and grp(surh)[1] is not None and grp(surl)[1] is not None:
        a, b = list(a)[::-1], list(b)[::-1]
    return vec(_clab.hzdsel(s, h, l, a, b))
