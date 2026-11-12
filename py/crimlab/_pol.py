# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Policing (R/pol.R)."""

from __future__ import annotations

from typing import Any

from . import _clab
from ._out import Record, Table, Vector, rec, tab, vec
from ._util import Vec, named


def clrsmp(clra: Vec, reca: Vec, clrb: Vec, recb: Vec) -> Record:
    """Clearance-rate decomposition and Simpson's paradox.

    Splits the change in the clearance rate between two areas or services
    (A and B) into a rate effect, from each offence type clearing more or
    less, and a mix effect, from the workload moving between types that
    clear at different rates (Kitagawa's decomposition); the two always add
    up to the change. Each part's mix term is measured from the mean of the
    two aggregate rates. See ``proof("clrsmp")``.

    Parameters
    ----------
    clra : offences cleared in A, by type, non-negative, optionally named.
    reca : offences recorded in A, by type, positive.
    clrb : offences cleared in B, by type, non-negative.
    recb : offences recorded in B, by type, positive.

    Returns
    -------
    ``parts``, a data frame with one row per type: ``part``, ``ratea``,
    ``rateb`` (clearance rates), ``sharea``, ``shareb`` (shares of the
    recorded workload), ``rateeff`` and ``mixeff``; and ``total``:
    ``ratea``, ``rateb`` (aggregate rates), ``change``, ``rateeff``,
    ``mixeff`` (``change = rateeff + mixeff``), ``up`` and ``down`` (types
    whose rate rose or fell) and ``simpson``: 2 when every type moved one
    way and the aggregate the other, 1 when the aggregate change and the
    rate effect differ in sign, 0 otherwise.

    Examples
    --------
    >>> r = clrsmp(clra={"minor": 90, "serious": 20}, reca=[100, 50],
    ...            clrb=[19, 60], recb=[20, 140])
    >>> r["parts"]["part"], round(r["total"]["change"], 4), r["total"]["simpson"]
    (['minor', 'serious'], -0.2396, 2.0)
    """
    (a, b, c, d), g = named("clrsmp", [("clra", clra), ("reca", reca), ("clrb", clrb),
                                       ("recb", recb)])
    return rec(_clab.clrsmp(a, b, c, d, g))


def expbnd(y: Vec, m: Vec, gamma: float, ref: Any = 1) -> Table:
    """Rates and disparities when exposure is measured by a proxy.

    If each group's true exposure lies within a factor ``gamma`` of an
    observed proxy, each group's true rate lies within a factor ``gamma``
    of its proxy rate, and the ratio of two groups' true rates within a
    factor ``gamma**2`` of the proxy ratio. The direction of a disparity is
    settled only when the proxy ratio lies outside that band. See
    ``proof("expbnd")``.

    Parameters
    ----------
    y : recorded counts by group, positive, optionally named.
    m : the exposure proxy by group, positive.
    gamma : the factor within which true exposure lies around the proxy,
        one number, at least 1.
    ref : the reference group, by name or by position from 1.

    Returns
    -------
    A data frame in the order of ``y``: ``group``, ``ref`` (True for the
    reference group), ``rate`` (``y / m``), ``lo`` and ``hi`` (the
    bounds on the true rate), ``ratio`` (the proxy rate against the
    reference group), ``rlo`` and ``rhi`` (the bounds on the true ratio)
    and ``sign``: 1 when the true ratio exceeds 1 for every admissible
    exposure, -1 when it is below 1 for every one, 0 when the data leave it
    open.

    Examples
    --------
    >>> r = expbnd(y={"white": 100, "black": 300}, m=[1000, 1000], gamma=1.5)
    >>> [round(x, 4) for x in r["ratio"]], r["sign"]
    ([1.0, 3.0], [0, 1])
    """
    (a, b), g = named("expbnd", [("y", y), ("m", m)])
    return tab(_clab.expbnd(a, b, gamma, ref, g))


def misrd(r1: float, r0: float, prev: float, se: float, sp: float) -> Vector:
    """Risk difference when group membership is misrecorded.

    When group membership is recorded with error at random, the measured
    risk difference is the true one times the factor
    ``(se + sp - 1) p (1 - p) / (d1 d0)``, never larger in size; with
    ``se + sp < 1`` it has the opposite sign. See ``proof("misrd")``.

    Parameters
    ----------
    r1, r0 : the true risks in groups 1 and 0, each from 0 to 1.
    prev : the true share in group 1, strictly between 0 and 1.
    se : sensitivity, from 0 to 1.
    sp : specificity, from 0 to 1.

    Returns
    -------
    ``true`` (``r1 - r0``), ``measured`` (the difference between the
    recorded groups), ``factor`` (``measured / true`` whenever ``true`` is
    not 0), ``d1`` and ``d0`` (the shares recorded in each group).

    Examples
    --------
    >>> {k: round(v, 4) for k, v in misrd(0.3, 0.1, 0.4, 0.8, 0.9).items()}
    {'true': 0.2, 'measured': 0.1426, 'factor': 0.7131, 'd1': 0.38, 'd0': 0.62}
    """
    return vec(_clab.misrd(r1, r0, prev, se, sp))


def rtmhs(x1: Vec, x2: Vec, cut: float, w: Vec | None = None) -> Vector:
    """Regression to the mean at selected hot spots.

    Places selected for a high count in one period tend to fall in the
    next even with no treatment. This reports the change at the selected
    places next to the same statistic on the data with its period-swapped
    copy, which is never positive; ``excess`` is the part of the fall that
    selection alone does not explain. See ``proof("rtmhs")``.

    Parameters
    ----------
    x1, x2 : counts by place in the first and second periods; when ``x1``,
        ``x2`` (and ``w``) are all named with unique names, they are matched
        by name.
    cut : the selection cut: places with ``x1 > cut`` are selected.
    w : None for equal weights, or one non-negative weight per place.

    Returns
    -------
    ``nsel`` (the weight selected), ``change`` (the weighted mean of
    ``x2 - x1`` at the selected places), ``mirror`` (the mean of ``x1 - x2``
    at places with ``x2 > cut``), ``sym`` (the change in the data with its
    swapped copy; never positive), ``excess`` (``change - sym``), ``low``
    (the mean change at places with ``x1 < cut``) and ``lowsym`` (its
    swapped-copy version; never negative). NaN where no place qualifies.

    Examples
    --------
    >>> {k: round(v, 4) for k, v in rtmhs([10, 2, 7], [6, 3, 8], 6).items()}
    {'nsel': 2.0, 'change': -1.5, 'mirror': -1.0, 'sym': -1.3333,
     'excess': -0.1667, 'low': 1.0, 'lowsym': 0.0}
    """
    (a, b, c), _ = named("rtmhs", [("x1", x1), ("x2", x2), ("w", w)], strict=False)
    return vec(_clab.rtmhs(a, b, cut, c))
