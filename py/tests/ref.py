# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Reading tests/data/ref.json, the values the R package gives
(dev/pyref.R), and comparing results with them."""

import json
import math
from pathlib import Path

CASES = json.loads((Path(__file__).parent / "data" / "ref.json").read_text(encoding="utf-8"))

# R's argument `or` is a Python keyword
RENAME: dict[str, str] = {}


def el(kind, e):
    if e is None:
        return None
    if kind == "num":
        return float(e)
    if kind == "int":
        return int(e)
    return e


def named(names, values):
    """Names and values as the package takes them: a dict, or a list of
    (name, value) tuples where a name repeats (a dict cannot)."""
    if len(set(names)) < len(names):
        return list(zip(names, values))
    return dict(zip(names, values))


def dec(x):
    """An R value as the Python argument that stands for it, in built-in
    types only: a named vector is a dict (a list of (name, value) tuples
    when a name repeats), a vector of one a scalar, a matrix a list of rows
    (a mapping from row name to row with row names; a mapping of columns,
    each a mapping from row name, with row and column names, both unique),
    a data frame a dict of columns (each named by the row names when it has
    them)."""
    if x is None:
        return None
    if "mat" in x:
        m = x["mat"]
        r, c = m["dim"]
        v = [float(e) if isinstance(e, (int, float)) else e for e in m["v"]]
        rows = [[v[i + j * r] for j in range(c)] for i in range(r)]
        rn, cn = m["rn"], m.get("cn")
        if rn and cn and len(set(rn)) == len(rn) and len(set(cn)) == len(cn):
            # row and column names: a mapping of columns, each by row name
            return {k: dict(zip(rn, (e[j] for e in rows))) for j, k in enumerate(cn)}
        return named(rn, rows) if rn else rows
    if "df" in x:
        cols = {k: dec(v) for k, v in x["df"].items()}
        if x["rn"]:
            cols = {k: named(x["rn"], v) for k, v in cols.items()}
        return cols
    kind = next(k for k in x if k != "names")
    v = [el(kind, e) for e in x[kind]]
    if x["names"] is not None:
        return named(x["names"], v)
    return v[0] if len(v) == 1 else v


def args(case):
    a = {RENAME.get(k, k): dec(v) for k, v in case["args"].items()}
    return a


class Tally:
    exact = 0
    close = 0


def near(a, b):
    if a == b:
        Tally.exact += 1
        return True
    Tally.close += 1
    return abs(a - b) <= 1e-12 * max(abs(a), abs(b))


def one(a, b, kind, where):
    if b is None:
        assert a is None, where
    elif kind == "num":
        b = float(b)
        assert isinstance(a, float), (where, a)
        if math.isnan(b):
            assert math.isnan(a), (where, a)
        elif math.isinf(b):
            assert a == b, (where, a, b)
        else:
            assert near(a, b), (where, a, b)
    elif kind == "int":
        assert type(a) is int and a == b, (where, a, b)
    elif kind == "lgl":
        assert type(a) is bool and a == b, (where, a, b)
    else:
        assert a == b, (where, a, b)


def vec(py, r, where):
    kind = next(k for k in r if k != "names")
    v = r[kind]
    if r["names"] is not None:
        assert isinstance(py, dict) and list(py) == r["names"], (where, py)
        py = list(py.values())
    elif not isinstance(py, list):
        py = [py]
    assert len(py) == len(v), (where, py, v)
    for i, (a, b) in enumerate(zip(py, v)):
        one(a, b, kind, f"{where}[{i}]")


def same(py, r, where=""):
    """py, a Python result, equals r, the R one: names, types, and every
    number to 1e-12 relative (Tally counts those equal exactly)."""
    for k in ("df", "list"):
        if k in r:
            assert isinstance(py, dict) and list(py) == list(r[k]), (where, list(py))
            for n, v in r[k].items():
                same(py[n], v, f"{where}${n}")
            return
    vec(py, r, where)
