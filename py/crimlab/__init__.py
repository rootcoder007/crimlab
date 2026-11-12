# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""crimlab: criminology and sociolegal statistics.

The same compiled kernels as the R package crimlab, with the same function
and argument names and the same results. Results are plain Python objects:
a data frame is a dict of columns (lists), a named vector a dict of
floats; :func:`frame` and :func:`arrays` turn them into pandas or NumPy
objects when those are installed.
"""

from __future__ import annotations

import importlib
from typing import Any

from . import _clab
from ._corr import hzdsel, incap, incapr, losbnd, sellab, sntbnd
from ._disp import dspbm, dyadrt, orrisk, outtst, vodtst
from ._frns import cllr, cllrmin, dbsrch, eceplt, gtmp, kinidx, lrcal, lrpost, mixlr, patidx, tipplt
from ._pol import clrsmp, expbnd, misrd, rtmhs
from ._proof import proof
from ._version import __version__

__all__ = [
    "__version__",
    "arrays",
    "cllr",
    "cllrmin",
    "clrsmp",
    "dbsrch",
    "dspbm",
    "dyadrt",
    "eceplt",
    "expbnd",
    "frame",
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


def frame(x: dict[str, Any]) -> Any:
    """A result as pandas: a dict of columns becomes a DataFrame, a dict of
    numbers a Series. Needs pandas.

    Examples
    --------
    >>> frame(orrisk(3, 0.4, 0.5)).index.tolist()  # doctest: +SKIP
    ['exposed', 'unexposed', 'rr', 'or']
    """
    pd = importlib.import_module("pandas")
    if all(isinstance(v, list) for v in x.values()):
        return pd.DataFrame(x)
    return pd.Series(x, dtype=float)


def arrays(x: Any) -> Any:
    """A result as NumPy: every list becomes an array (None, R's NA, is
    NaN in a numeric one), nested dicts alike. Needs NumPy.

    Examples
    --------
    >>> arrays(vodtst([6, 4], [5, 5]))["ratio"]  # doctest: +SKIP
    array([1.        , 0.66666667])
    """
    np = importlib.import_module("numpy")
    if isinstance(x, dict):
        return {k: arrays(v) for k, v in x.items()}
    if isinstance(x, list):
        num = all(v is None or isinstance(v, (int, float)) for v in x)
        return np.array(x, dtype=float) if num and None in x else np.asarray(x)
    return x
