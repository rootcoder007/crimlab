# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Parity with R: every case of tests/data/ref.json (written by
dev/pyref.R) gives the R result, or stops with R's message."""

import pytest
from ref import CASES, Tally, args, needs, same

import crimlab

FNS = sorted({c["fn"] for c in CASES})


def test_every_exported_function_has_cases():
    extra = {"selftest", "frame", "arrays", "__version__"}
    assert set(FNS) == set(crimlab.__all__) - extra
    assert len(FNS) == 27


@pytest.mark.parametrize("case", [c for c in CASES if "out" in c],
                         ids=lambda c: c["fn"])
def test_result_matches_r(case):
    if needs(case):
        pytest.importorskip(needs(case))
    same(getattr(crimlab, case["fn"])(**args(case)), case["out"], case["fn"])


@pytest.mark.parametrize("case", [c for c in CASES if "err" in c],
                         ids=lambda c: c["fn"])
def test_error_matches_r(case):
    if case["fn"] == "proof" and case["args"]["fn"] == {"chr": [None], "names": None}:
        pytest.skip("proof(NA): None is proof()'s default in Python; test_misc checks [None]")
    if needs(case):
        pytest.importorskip(needs(case))
    with pytest.raises((ValueError, TypeError)) as e:
        getattr(crimlab, case["fn"])(**args(case))
    assert str(e.value) == case["err"]


def test_tally():
    # every number above was compared; most are bit-for-bit R's
    print(f"\nnumbers equal to R exactly: {Tally.exact}, within 1e-12: {Tally.close}")
    assert Tally.exact > 0
