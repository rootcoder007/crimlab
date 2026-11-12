# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""crimlab: criminology and sociolegal statistics.

The same compiled kernels as the R package crimlab, with the same function
and argument names and the same results. Results are plain Python objects:
a data frame is a Table, a dict of columns (lists); a named vector a
Vector, a dict of floats; a list a Record. They print as R prints them and
export themselves as CSV or JSON; a Table also sorts, filters and joins.
The package imports nothing outside the standard library.
"""

from __future__ import annotations

from typing import Any

from . import _clab
from ._corr import hzdsel, incap, incapr, losbnd, sellab, sntbnd
from ._disp import dspbm, dyadrt, orrisk, outtst, vodtst
from ._frns import cllr, cllrmin, dbsrch, eceplt, gtmp, kinidx, lrcal, lrpost, mixlr, patidx, tipplt
from ._out import Record, Table, Vector
from ._pol import clrsmp, expbnd, misrd, rtmhs
from ._proof import proof
from ._version import __version__

__all__ = [
    "Record",
    "Table",
    "Vector",
    "__version__",
    "cllr",
    "cllrmin",
    "clrsmp",
    "dbsrch",
    "dspbm",
    "dyadrt",
    "eceplt",
    "expbnd",
    "gtmp",
    "hzdsel",
    "incap",
    "incapr",
    "kinidx",
    "losbnd",
    "lrcal",
    "lrpost",
    "misrd",
    "mixlr",
    "orrisk",
    "outtst",
    "patidx",
    "proof",
    "rtmhs",
    "selftest",
    "sellab",
    "sntbnd",
    "tipplt",
    "vodtst",
]


def selftest() -> dict[str, Any]:
    """Run the compiled self-test of the kernels (src/selftest.cpp).

    Returns
    -------
    ``n`` (checks run), ``nfail`` (checks failed) and ``report`` (one line
    per failed check; empty when all pass).

    Examples
    --------
    >>> r = selftest()
    >>> r["nfail"], r["report"]
    (0, '')
    """
    return _clab.selftest()
