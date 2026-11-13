# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The self-test, the proofs, the version, the result helpers, the plots,
the shapes and types of results, and the length hook."""

import gc
import math
import re
import subprocess
import sys
import tracemalloc
from pathlib import Path

import pytest
from ref import CASES, args


def rnames():
    """The functions the R package exports, or None outside the repository."""
    ns = Path(__file__).resolve().parents[2] / "NAMESPACE"
    if not ns.is_file():  # an sdist
        return None
    return set(re.findall(r"^export\((\w+)\)", ns.read_text(), re.MULTILINE))


import crimlab
from crimlab import _clab


def test_the_compiled_self_test_passes_every_check():
    r = crimlab.selftest()
    assert r["report"] == ""
    assert r["nfail"] == 0
    assert r["n"] >= 13700


def test_the_version_is_the_r_version():
    assert re.fullmatch(r"\d+(\.\d+)+", crimlab.__version__)
    desc = Path(__file__).resolve().parents[2] / "DESCRIPTION"
    if desc.is_file():  # in the repository, not in an sdist
        v = re.search(r"^Version:\s*(\S+)", desc.read_text(), re.MULTILINE).group(1)
        assert crimlab.__version__ == v.replace("-", ".")


def test_proof_names_real_theorems():
    p = crimlab.proof()
    assert list(p) == ["fn", "theorem", "file"]
    assert len(set(p["theorem"])) == len(p["theorem"])
    # proof() lists the theorems of every function of the R package, those
    # not bound in Python too: every name is an export of the R package
    r = rnames()
    if r is not None:
        assert set(p["fn"]) <= r
    assert len(crimlab.proof("dspbm")["fn"]) == 2
    root = Path(crimlab.__file__).parent
    for th, f in zip(p["theorem"], p["file"]):
        src = (root / f).read_text(encoding="utf-8").splitlines()
        nm = th.rsplit(".", 1)[1]
        assert any(re.match(rf"theorem {nm}\b", ln) for ln in src), nm
        assert f"namespace {th.rsplit('.', 1)[0]}" in src
    for f in root.glob("lean/*.lean"):
        assert "sorry" not in f.read_text(encoding="utf-8")
    two = crimlab.proof(["cllr", "dbsrch"])
    assert list(dict.fromkeys(two["fn"])) == ["cllr", "dbsrch"]
    with pytest.raises(TypeError, match="function names"):
        crimlab.proof([None])
    with pytest.raises(TypeError, match="function names"):
        crimlab.proof(3)


def test_lengths_past_int_max_are_refused():
    assert _clab.tlen(5.0, "x") == 5
    assert _clab.tlen(2147483647.0, False) == 2147483647
    with pytest.raises(ValueError, match=r"^clab_tlen: `x` is too long for crimlab: "
                                         r"more than 2147483647 values\.$"):
        _clab.tlen(2147483648.0, "x")
    with pytest.raises(ValueError, match="^clab_tlen: a size or an index is out of range$"):
        _clab.tlen(3e9, False)
    assert _clab.tlen(7.0, None) == "7"
    assert _clab.tlen(3e9, None) == "3e+09"
    with pytest.raises(TypeError):
        _clab.tlen("x", None)
    with pytest.raises(TypeError, match="expected 2 arguments, got 1"):
        _clab.tlen(1.0)


@pytest.mark.parametrize("x, s", [
    (1e5, "1e+05"), (123456, "123456"), (0.1 + 0.2, "0.3"), (1e-4, "1e-04"),
    (1 / 3, "0.333333333333333"), (2, "2"), (10, "10"), (1e15, "1e+15"), (1e-15, "1e-15"),
    (123456.7, "123456.7"), (5e6, "5e+06"), (100, "100"), (1234567890123456, "1234567890123456"),
    (0.000123, "0.000123"), (-2.5, "-2.5"), (1e100, "1e+100"), (1e-300, "1e-300"),
    (2.0**53, "9007199254740992"), (0, "0"), (math.inf, "Inf"), (-math.inf, "-Inf"),
    (math.nan, "NaN"),
])
def test_numbers_are_written_as_r_writes_them(x, s):
    # R 4.6: as.character(c(1e5, 123456, 0.1 + 0.2, ...))
    assert _clab.rchar(float(x)) == s


def test_results_have_the_documented_shapes_and_types():
    d = {"a": 600, "b": 400}
    n = {"a": 700, "b": 300}

    def shape(x):
        return " ".join(f"{k}:{type(v[0]).__name__}" for k, v in x.items())

    assert shape(crimlab.vodtst(d, n)) == \
        "group:str ref:bool shday:float shnight:float ratio:float sign:int"
    assert shape(crimlab.gtmp([0.1, 0.2], [0.3, None], [False, True])["locus"]) == \
        "locus:str hom:bool mp:float"
    assert shape(crimlab.expbnd(d, n, 1.5)) == \
        "group:str ref:bool rate:float lo:float hi:float ratio:float rlo:float rhi:float sign:int"
    r = crimlab.incap({"low": 2, "high": 10}, 0.1, 1, share=[0.8, 0.2])
    assert r["group"][-1] == "all" and r["marg"][-1] is None and r["q"][-1] is None
    s = crimlab.sntbnd([1, 0, 1, 0], ["x", "x", "y", "y"], "y")
    assert s["share"][2] is None
    p = crimlab.patidx([1, 2], [1, 2], [0.1, 0.2], mother=[3, 3])
    assert math.isnan(p["locus"]["idx"][0])


def test_no_third_party_package_is_imported():
    # in a fresh interpreter: import crimlab, run every reference case
    code = (
        "import sys, crimlab\n"
        "from ref import CASES, args\n"
        "for c in CASES:\n"
        "    try:\n"
        "        getattr(crimlab, c['fn'])(**args(c))\n"
        "    except (ValueError, TypeError):\n"
        "        pass\n"
        "crimlab.selftest(); crimlab.proof()\n"
        "print(sorted(m for m in ('pandas', 'numpy', 'matplotlib') if m in sys.modules))\n"
    )
    r = subprocess.run([sys.executable, "-c", code], cwd=Path(__file__).parent,
                       capture_output=True, text=True, check=False, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.strip() == "[]"


def test_no_memory_is_kept_by_calls_or_errors():
    # every reference case, valid and invalid, run many times: the memory
    # Python allocates must not grow with the number of runs
    def run(k):
        for _ in range(k):
            for c in CASES:
                try:
                    getattr(crimlab, c["fn"])(**args(c))
                except (ValueError, TypeError):
                    pass

    run(3)
    gc.collect()
    tracemalloc.start()
    run(20)
    gc.collect()
    a = tracemalloc.get_traced_memory()[0]
    run(60)
    gc.collect()
    b = tracemalloc.get_traced_memory()[0]
    tracemalloc.stop()
    assert b - a < 64 * 1024, (a, b)
