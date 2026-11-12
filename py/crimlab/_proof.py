# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The machine-checked proofs behind each function (R/proof.R)."""

from __future__ import annotations

from typing import Any

from ._thm import THM
from ._util import flat, isseq


def proof(fn: Any = None) -> dict[str, list[str]]:
    """Machine-checked proofs behind the package.

    Names the Lean theorem each function rests on. The theorems ship with
    the package: ``file`` is the path below the package's directory, and
    each file compiles against Mathlib.

    Parameters
    ----------
    fn : None for every function, or function names (one name or a list).

    Returns
    -------
    A data frame: ``fn``, ``theorem`` and ``file``.

    Examples
    --------
    >>> proof("vodtst")["theorem"]
    ['Crimlab.Veil.odds_ratio_identifies_bias', 'Crimlab.Veil.share_sign', 'Crimlab.Veil.exposure_shift_confounds']
    """
    rows = THM
    if fn is not None:
        fn = [fn] if isinstance(fn, str) else flat(fn)
        if not isseq(fn) or not all(isinstance(f, str) for f in fn):
            raise TypeError("proof: `fn` must be function names.")
        known = {t[0] for t in THM}
        miss = list(dict.fromkeys(f for f in fn if f not in known))
        if miss:
            raise ValueError(f"proof: no proof recorded for: {', '.join(miss)}.")
        rows = tuple(t for t in THM if t[0] in set(fn))
    return {
        "fn": [t[0] for t in rows],
        "theorem": [f"Crimlab.{t[1]}.{t[2]}" for t in rows],
        "file": [f"lean/{t[1]}.lean" for t in rows],
    }
