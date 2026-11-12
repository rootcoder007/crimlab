# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Results as tables of our own: printing, rows, CSV, JSON, and sort,
filter and join, each against a plain-Python reference."""

import csv
import io
import json
import math
import random
import time

import pytest

import crimlab
from crimlab import Record, Table, Vector


def test_results_are_still_dicts():
    r = crimlab.vodtst({"w": 600, "b": 400}, {"w": 700, "b": 300})
    assert isinstance(r, Table) and isinstance(r, dict)
    assert r["ratio"][0] == 1.0 and list(r) == ["group", "ref", "shday", "shnight", "ratio", "sign"]
    assert isinstance(crimlab.orrisk(3, 0.4, 0.5), Vector)
    g = crimlab.gtmp([0.1, 0.2], [0.3, None], [False, True])
    assert isinstance(g, Record) and isinstance(g["locus"], Table)
    assert isinstance(crimlab.clrsmp([9, 2], [10, 5], [1, 6], [2, 14])["total"], Vector)
    assert isinstance(crimlab.proof("cllr"), Table)
    t = Table({"sort": [1, 2], "keys": [3, 4]})
    assert t["sort"] == [1, 2] and list(t.keys()) == ["sort", "keys"]
    assert t.get("x") is None and dict(t.items()) == {"sort": [1, 2], "keys": [3, 4]}


def test_printing():
    t = Table({"g": ["a", "bb"], "x": [1.5, None], "y": [float("nan"), float("inf")],
               "b": [True, False], "n": [1, 22]})
    assert repr(t) == ("   g   x   y     b  n\n"
                       "1  a 1.5 NaN  TRUE  1\n"
                       "2 bb  NA Inf FALSE 22")
    assert str(Table({"x": []})) == " x\n(0 rows)"
    long = Table({"i": list(range(100))})
    s = repr(long)
    assert s.splitlines()[0] == "      i" or s.splitlines()[0].strip() == "i"
    assert "..." in s and s.endswith("(100 rows)") and len(s.splitlines()) == 1 + 10 + 1 + 5 + 1
    assert repr(Vector({"a": 0.5, "bb": 0.25})) == "   a   bb\n0.50 0.25"
    assert repr(Table({"x": [1e-10, 1.0]})).splitlines()[1].endswith("1e-10")
    assert "$log10mp\n[1] " in repr(crimlab.gtmp([0.1], [0.2]))


def test_rows_row_and_shape():
    t = Table({"a": [1, 2, 3], "b": ["x", "y", "z"]})
    assert t.shape == (3, 2) and Table().shape == (0, 0)
    assert t.rows() == [{"a": 1, "b": "x"}, {"a": 2, "b": "y"}, {"a": 3, "b": "z"}]
    assert t.row(-1) == {"a": 3, "b": "z"}
    with pytest.raises(IndexError):
        t.row(3)
    with pytest.raises(ValueError, match="one length"):
        Table({"a": [1], "b": [1, 2]})


def test_csv_reads_back():
    t = Table({"g": ["a,b", 'q"t', "line\nbreak"], "x": [0.1 + 0.2, None, float("-inf")],
               "b": [True, False, True]})
    s = t.csv()
    assert s.startswith("g,x,b\r\n") and s.endswith("\r\n")
    back = list(csv.reader(io.StringIO(s, newline="")))
    assert back[1:] == [["a,b", "0.30000000000000004", "TRUE"], ['q"t', "", "FALSE"],
                        ["line\nbreak", "-Inf", "TRUE"]]
    assert float(back[1][1]) == 0.1 + 0.2


def test_json_reads_back():
    r = crimlab.incap({"low": 2, "high": 10}, 0.1, 1, share=[0.8, 0.2])
    d = json.loads(r.json())
    assert d["group"] == ["low", "high", "all"] and d["marg"][2] is None
    assert d["lam"] == r["lam"]
    assert json.loads(Table({"x": [float("nan"), math.inf, -math.inf]}).json()) == \
        {"x": ["NaN", "Inf", "-Inf"]}
    g = json.loads(crimlab.gtmp([0.1], [0.2]).json())
    assert set(g) == {"locus", "log10mp", "log10lr"}
    assert json.loads(crimlab.orrisk(3, 0.4, 0.5).json())["or"] == 3.0


# References in plain Python

def miss(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def refsort(t, by, desc):
    n = len(next(iter(t.values())))
    idx = list(range(n))
    for c, d in reversed(list(zip(by, desc))):  # stable: last key first
        have = [i for i in idx if not miss(t[c][i])]
        none = [i for i in idx if miss(t[c][i])]
        have.sort(key=lambda i: t[c][i], reverse=d)
        idx = have + none
    return {k: [v[i] for i in idx] for k, v in t.items()}


def rnd(n, seed):
    g = random.Random(seed)
    return Table({
        "k": [g.choice(["a", "b", "c", None]) for _ in range(n)],
        "x": [g.choice([1.0, 2.0, 2.0, 3.5, None, math.nan]) for _ in range(n)],
        "i": list(range(n)),
    })


def same(a, b):
    def norm(v):
        return "NaN" if isinstance(v, float) and math.isnan(v) else v

    return {k: [norm(v) for v in c] for k, c in a.items()} == \
        {k: [norm(v) for v in c] for k, c in b.items()}


@pytest.mark.parametrize("seed", range(12))
def test_sort_against_a_reference(seed):
    t = rnd(random.Random(seed).randint(0, 40), seed)
    before = {k: list(v) for k, v in t.items()}
    for by, desc in ((["x"], [False]), (["x"], [True]), (["k", "x"], [False, True]),
                     (["k", "x"], [True, False]), (["k"], [True])):
        assert same(t.sort(by, desc), refsort(t, by, desc)), (by, desc)
    assert t == before or same(t, before)
    # stable: equal keys keep their order (i ascending)
    s = t.sort("k")
    for a, b in zip(s.rows(), s.rows()[1:]):
        if a["k"] == b["k"]:
            assert a["i"] < b["i"]


def test_sort_errors_and_text_order():
    t = Table({"s": ["b", "B", "a", None], "m": [1, "x", 2, 3]})
    assert t.sort("s")["s"] == ["B", "a", "b", None]
    assert t.sort("s", desc=True)["s"] == ["b", "a", "B", None]
    with pytest.raises(TypeError, match="sort: column `m` mixes text and numbers"):
        t.sort("m")
    with pytest.raises(ValueError, match="sort: no column `z`"):
        t.sort("z")
    with pytest.raises(ValueError, match="one per column"):
        t.sort(["s"], desc=[True, False])
    assert Table({"x": []}).sort("x") == {"x": []}
    assert Table({"x": [1]}).sort("x") == {"x": [1]}


OPS = {"==": lambda a, b: a == b, "!=": lambda a, b: a != b, "<": lambda a, b: a < b,
       "<=": lambda a, b: a <= b, ">": lambda a, b: a > b, ">=": lambda a, b: a >= b,
       "in": lambda a, b: a in b}


@pytest.mark.parametrize("seed", range(8))
def test_filter_against_a_reference(seed):
    t = rnd(random.Random(seed).randint(0, 40), seed)
    for op, val in (("==", 2.0), ("!=", 2.0), ("<", 2.0), ("<=", 2.0), (">", 1.0),
                    (">=", 3.5), ("in", [1.0, 3.5])):
        keep = [not miss(v) and OPS[op](v, val) for v in t["x"]]
        want = {k: [a for a, b in zip(c, keep) if b] for k, c in t.items()}
        assert same(t.filter(x=(op, val)), want), op
    keep = [k == "a" and not miss(x) and x > 1 for k, x in zip(t["k"], t["x"])]
    want = {k: [a for a, b in zip(c, keep) if b] for k, c in t.items()}
    assert same(t.filter(k="a", x=(">", 1)), want)
    assert same(t.filter(lambda r: r["i"] % 2 == 0, k=("in", ["a", "b"])),
                {k: [a for a, r in zip(c, t.rows()) if r["i"] % 2 == 0 and r["k"] in ("a", "b")]
                 for k, c in t.items()})


def test_filter_errors():
    t = Table({"g": ["a", "b"], "x": [1.0, None], "ref": [True, False]})
    assert t.filter(ref=False)["g"] == ["b"]
    assert t.filter(x=("!=", 5))["g"] == ["a"]  # a missing value never passes
    with pytest.raises(ValueError, match="filter: no column `z`"):
        t.filter(z=1)
    with pytest.raises(ValueError, match="filter: unknown operator `~`"):
        t.filter(x=("~", 1))
    with pytest.raises(TypeError, match="filter: column `g` cannot be compared"):
        t.filter(g=("<", 1))
    with pytest.raises(TypeError, match="`in` for column `g` needs a list"):
        t.filter(g=("in", "a"))
    with pytest.raises(TypeError, match="function of a row"):
        t.filter(3)


def refjoin(a, b, on, how):
    ka = [tuple(a[c][i] for c in on) for i in range(len(a[on[0]]))]
    kb = [tuple(b[c][i] for c in on) for i in range(len(b[on[0]]))]

    def ok(k):
        return not any(miss(v) for v in k)

    pairs, used = [], set()
    for i, k in enumerate(ka):
        m = [j for j, kk in enumerate(kb) if ok(k) and kk == k]
        pairs += [(i, j) for j in m]
        used.update(m)
        if not m and how in ("left", "full"):
            pairs.append((i, -1))
    if how in ("right", "full"):
        pairs += [(-1, j) for j in range(len(kb)) if j not in used]
    out = {c: [a[c][i] if i >= 0 else b[c][j] for i, j in pairs] for c in on}
    for c in a:
        if c not in on:
            out[c + ".x" if c in b else c] = [a[c][i] if i >= 0 else None for i, _ in pairs]
    for c in b:
        if c not in on:
            out[c + ".y" if c in a else c] = [b[c][j] if j >= 0 else None for _, j in pairs]
    return out


@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("how", ["inner", "left", "right", "full"])
def test_join_against_a_reference(seed, how):
    g = random.Random(seed)
    a = rnd(g.randint(0, 30), seed)
    b = rnd(g.randint(0, 30), seed + 100)
    for on in (["k"], ["k", "x"]):
        assert same(a.join(b, on, how), refjoin(a, b, on, how)), (on, how)
    assert same(a.join(dict(b), "k", how), refjoin(a, b, ["k"], how))


def test_join_rules():
    a = Table({"k": ["a", "a", "b", None], "v": [1, 2, 3, 4]})
    b = {"k": ["a", "a", "c", None], "v": [10, 20, 30, 40]}
    r = a.join(b, "k")
    assert r.shape == (4, 3) and list(r) == ["k", "v.x", "v.y"]  # 2 x 2 pairings
    assert a.join(b, "k", "full").shape[0] == 4 + 2 + 2  # pairs, b and None left, c and None right
    # numbers match numbers (1 == 1.0), never text or bools
    assert Table({"k": [1]}).join({"k": [1.0]}, "k").shape[0] == 1
    assert Table({"k": [1]}).join({"k": ["1"]}, "k").shape[0] == 0
    assert Table({"k": [1]}).join({"k": [True]}, "k").shape[0] == 0
    with pytest.raises(ValueError, match="no key column `z` in this table"):
        a.join({"k": [1], "q": [2]}, ["k", "z"])
    with pytest.raises(ValueError, match="`how` must be inner, left, right or full"):
        a.join(b, "k", "outer")
    with pytest.raises(ValueError, match="one length"):
        a.join({"k": [1], "q": [1, 2]}, "k")


def test_join_is_linear():
    def run(n):
        a = Table({"k": list(range(n)), "x": [1.0] * n})
        b = Table({"k": list(range(n - 1, -1, -1)), "y": [2.0] * n})
        t0 = time.perf_counter()
        r = a.join(b, "k")
        assert r.shape == (n, 3)
        return time.perf_counter() - t0

    run(1000)
    small, big = min(run(20000) for _ in range(3)), min(run(100000) for _ in range(3))
    # five times the rows: about five times the time; quadratic would be 25
    assert big < 12 * small, (small, big)


def test_values_of_no_kind_and_control_characters():
    assert json.loads(Table({"s": ["a\x01b\t\n\\\""]}).json()) == {"s": ["a\x01b\t\n\\\""]}
    t = Table({"o": [object(), object()], "k": [1, 2]})
    assert "<object object" in repr(t) and t.csv().count("<object object") == 2
    with pytest.raises(TypeError, match="sort: column `o` holds values that are neither"):
        t.sort("o")
    with pytest.raises(TypeError, match="join: a key holds values that are neither"):
        t.join(t, "o")
