# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Parity with R through clab_jrun(): every case of cli/test.tsv (the
values of the R package, dev/clref.R) gives R's result in Python, or stops
with R's message. Outside the repository (an sdist) there is no table."""

import json
import keyword
import math
from pathlib import Path

import pytest

import crimlab
from crimlab import _j

# the functions of crimlab/_j.py
NJ = 166

TSV = Path(__file__).resolve().parents[2] / "cli" / "test.tsv"
# cli/tol.tsv: the functions with bounds other than 1e-12 relative
# (relative, absolute: a number within either passes)
TOL = {}
if TSV.with_name("tol.tsv").is_file():
    for ln in TSV.with_name("tol.tsv").read_text(encoding="utf-8").splitlines():
        if ln and not ln.startswith("#"):
            p = ln.split("\t")
            TOL[p[0]] = (float(p[1]), float(p[2]))
CASES = []
if TSV.is_file():
    for ln in TSV.read_text(encoding="utf-8").splitlines():
        if ln and not ln.startswith("#"):
            f, a, w = ln.split("\t")
            # the functions with their own entry points have tests/ref.py
            if hasattr(_j, f):
                CASES.append((f, a, w))


def nodup(pairs):
    keys = [k for k, _ in pairs]
    if len(set(keys)) < len(keys):
        pytest.skip("a name given twice: a dict cannot hold it (R and JSON can)")
    return dict(pairs)


def kw(a):
    json.loads(a, object_pairs_hook=nodup)
    return {(k + "_" if keyword.iskeyword(k) else k): v for k, v in json.loads(a).items()}


def num(x):
    return float(x) if isinstance(x, str) and x in ("NaN", "Inf", "-Inf") else x


def same(py, r, where, rel=1e-12, ab=0.0):
    r = num(r)
    if isinstance(r, dict):
        assert isinstance(py, dict) and list(py) == list(r), (where, py)
        for k in r:
            same(py[k], r[k], f"{where}${k}", rel, ab)
    elif isinstance(r, list):
        assert isinstance(py, list) and len(py) == len(r), (where, py)
        for i, (a, b) in enumerate(zip(py, r)):
            same(a, b, f"{where}[{i}]", rel, ab)
    elif isinstance(r, bool) or r is None or isinstance(r, str):
        assert py == r and type(py) is type(r), (where, py, r)
    elif math.isnan(r):
        assert isinstance(py, float) and math.isnan(py), (where, py)
    else:
        assert isinstance(py, (int, float)) and not isinstance(py, bool), (where, py)
        tol = rel if r == 0 else rel * max(abs(py), abs(r))
        assert py == r or abs(py - r) <= max(tol, ab), (where, py, r)


@pytest.mark.parametrize("f,a,w", [c for c in CASES if not c[2].startswith("!")], ids=lambda x: x[:20])
def test_result_matches_r(f, a, w):
    same(getattr(crimlab, f)(**kw(a)), json.loads(w), f, *TOL.get(f, (1e-12, 0.0)))


@pytest.mark.parametrize("f,a,w", [c for c in CASES if c[2].startswith("!")], ids=lambda x: x[:20])
def test_error_matches_r(f, a, w):
    with pytest.raises((ValueError, TypeError)) as e:
        getattr(crimlab, f)(**kw(a))
    assert str(e.value) == w[1:]


def test_every_function_of_the_package_has_cases():
    if not CASES:
        pytest.skip("no cli/test.tsv outside the repository")
    fns = {f for f in dir(_j) if callable(getattr(_j, f)) and getattr(_j, f).__module__ == _j.__name__}
    assert fns == {c[0] for c in CASES}
    assert fns <= set(crimlab.__all__)
    assert len(fns) == NJ
