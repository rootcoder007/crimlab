# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Reading tests/data/ref.json, the values the R package gives
(dev/pyref.R), and comparing results with them."""

import json
import math
from pathlib import Path

import pytest

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


def dec(x):
    """An R value as the Python argument that stands for it: a named vector
    is a dict, a vector of one a scalar, a matrix a list of rows (a dict of
    rows when it has row names; a pandas DataFrame when it has row and
    column names), a data frame a dict of columns (each a dict by row name
    when it has row names)."""
    if x is None:
        return None
    if "mat" in x:
        m = x["mat"]
        r, c = m["dim"]
        v = [float(e) if isinstance(e, (int, float)) else e for e in m["v"]]
        rows = [[v[i + j * r] for j in range(c)] for i in range(r)]
        if m["rn"] and m.get("cn"):
            pd = pytest.importorskip("pandas")
            return pd.DataFrame(rows, index=m["rn"], columns=m["cn"])
        return dict(zip(m["rn"], rows)) if m["rn"] else rows
    if "df" in x:
        cols = {k: dec(v) for k, v in x["df"].items()}
        if x["rn"]:
            cols = {k: dict(zip(x["rn"], v)) for k, v in cols.items()}
        return cols
    kind = next(k for k in x if k != "names")
    v = [el(kind, e) for e in x[kind]]
    nm = x["names"]
    if nm is not None and len(set(nm)) < len(nm):  # a dict cannot repeat a name
        pd = pytest.importorskip("pandas")
        return pd.Series(v, index=nm)
    if nm is not None:
        return dict(zip(nm, v))
    return v[0] if len(v) == 1 else v


def needs(case):
    """The optional package a case needs to be written in Python: pandas
    for a table with row and column names or a name given twice."""
    def dup(x):
        if isinstance(x, dict):
            nm = x.get("names")
            if isinstance(nm, list) and len(set(nm)) < len(nm):
                return True
            return any(dup(v) for v in x.values())
        return False

    if '"cn": [' in json.dumps(case["args"]) or dup(case["args"]):
        return "pandas"
    return None


def args(case):
    a = {RENAME.get(k, k): dec(v) for k, v in case["args"].items()}
    if case["fn"] in ("tipplt", "eceplt"):
        a["draw"] = False
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
