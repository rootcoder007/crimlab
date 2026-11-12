# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The edge cases of the R test suite (tests/testthat) that are not plain
calls with a reference value: identities, exact values, limits and the
checks behind the messages."""

import math
import warnings

import pytest

from crimlab import (
    cllr,
    cllrmin,
    dbsrch,
    dyadrt,
    eceplt,
    expbnd,
    gtmp,
    hzdsel,
    incap,
    incapr,
    kinidx,
    losbnd,
    lrcal,
    lrpost,
    mixlr,
    orrisk,
    patidx,
    rtmhs,
    sntbnd,
    tipplt,
    vodtst,
)


def test_vodtst_and_dyadrt_label_what_they_return():
    assert "adj" not in vodtst({"w": 60, "b": 40}, {"w": 50, "b": 50})
    pop = {"a": 800, "b": 150, "c": 50}
    count = [[pop[i] * pop[j] * 0.002 for j in pop] for i in pop]
    r = dyadrt(count, pop)
    assert r["off"][1] == "b" and r["vic"][1] == "a"


def test_relative_risk_between_1_and_the_odds_ratio():
    for o in (0.05, 0.5, 1, 3, 400):
        for q in (0.01, 0.4, 0.99):
            for s in (0.02, 0.5, 0.9):
                r = orrisk(o, q, s)
                assert abs(math.log(r["rr"])) <= abs(math.log(o)) + 1e-12


def test_gtmp_ratio_is_the_inverse_and_loci_are_labelled():
    r = gtmp({"a": 0.1, "b": 0.2}, [0.3, None], hom=[False, True])
    assert r["log10lr"] == -r["log10mp"]
    assert r["locus"]["locus"] == ["a", "b"]
    p, q = 0.07, 0.11
    assert gtmp(p, q)["log10lr"] == pytest.approx(-math.log10(2 * p * q), rel=1e-10)
    assert gtmp(p, hom=True)["log10lr"] == pytest.approx(-math.log10(p * p), rel=1e-10)


def test_q_is_missing_in_every_form():
    ref = gtmp([0.2, 0.3], hom=True)
    for q in ([None, None], [math.nan, None], (None, math.nan)):
        assert gtmp([0.2, 0.3], q, hom=True) == ref
    assert gtmp(0.2, 0, hom=True) == gtmp(0.2, hom=True)
    with pytest.raises(TypeError, match="`q` must be numbers"):
        gtmp([0.2, 0.3], [None, "a"], hom=True)
    with pytest.raises(ValueError, match="`q`"):
        gtmp(0.2, 1)


def test_lrpost_applies_bayes_rule():
    assert lrpost(1e6, 1 / 5e6)[0] == pytest.approx((1e6 / (5e6 - 1)) / (1 + 1e6 / (5e6 - 1)),
                                                     rel=1e-10)
    assert lrpost(1, 0.3) == [0.3]
    assert lrpost([1, 3], [0.5, 0.25]) == [0.5, 0.5]
    assert lrpost(1e308, 0.5)[0] == pytest.approx(1, rel=1e-10)
    assert lrpost(1e300, 1e-300) == [0.5]


def test_kinidx_sums_the_loci():
    p, q = 0.1, 0.2
    g1 = {"D3": [1, 2], "vWA": [1, 1]}
    r = kinidx(g1, g1, [[p, q], [p, p]])
    assert r["locus"]["locus"] == ["D3", "vWA"]
    want = math.log10((p + q) / (4 * p * q)) + math.log10(1 / p)
    assert r["log10idx"] == pytest.approx(want, rel=1e-10)
    assert r["post"] == pytest.approx(1 / (1 + 10 ** -r["log10idx"]), rel=1e-10)


def test_mixlr_forms_the_ratio_of_its_two_sides():
    r = mixlr({"a": 0.1, "b": 0.2, "c": 0.05, "d": 0.15}, [], 0, ["c", "d"], 1)
    assert list(r.values()) == pytest.approx([1, 2 * 0.05 * 0.15, 1 / (2 * 0.05 * 0.15)],
                                             rel=1e-10)
    assert mixlr({"a": 0.1, None: 0.2}, "a", 1, "a", 1)["lr"] == 1


def test_dbsrch_bounds_the_exact_chance_by_np():
    r = dbsrch(1e-3, 100)
    assert r["np"] == 0.1
    assert r["any"] == pytest.approx(1 - (1 - 1e-3) ** 100, rel=1e-10)
    assert r["any"] <= r["np"]
    assert dbsrch(1e-6, 1e6, pop=1e6)["post"] == 1
    assert dbsrch(1e-20, 10)["any"] == pytest.approx(1e-19, rel=1e-10)


def test_calibrated_ratios_are_exact():
    s, d = [3, 0.7, 20, 5], [0.2, 1, 4, 0.05]
    cal = lrcal(s, d)
    assert cal == {"same": [1.0, 1.0, math.inf, math.inf], "diff": [0.0, 1.0, 1.0, 0.0]}
    a = cllrmin(s, d)

    def fin(x):
        return [min(max(v, 1e-300), 1e300) for v in x]

    assert cllr(fin(cal["same"]), fin(cal["diff"])) == pytest.approx(a["min"], rel=1e-10)
    full = 0.5 * (sum(math.log2(1 + 1 / v) for v in s) / 4 + sum(math.log2(1 + v) for v in d) / 4)
    assert a["cllr"] == pytest.approx(full, rel=1e-12)
    assert a["cal"] == pytest.approx(full - a["min"], rel=1e-12) and a["cal"] > 0
    r = cllrmin([50, 200, 8], [0.01, 0.2, 1.5])
    assert r["cal"] == r["cllr"]


def test_the_plots_return_what_they_draw():
    s, d = [50, 200, 8, 0.5], [0.01, 0.2, 1.5, 0.001]
    t = tipplt(s, d, draw=False)
    assert t["same"]["share"] == [1, 0.75, 0.5, 0.25]
    assert t["same"]["log10lr"] == pytest.approx(sorted(math.log10(v) for v in s), rel=1e-12)
    e = eceplt(s, d, lo=[-1, 0, 1], draw=False)
    m = cllrmin(s, d)
    assert e["ece"][1] == pytest.approx(m["cllr"], rel=1e-10)
    assert e["min"][1] == pytest.approx(m["min"], rel=1e-10)
    assert e["ref"][1] == pytest.approx(1, rel=1e-10)
    assert all(a <= b + 1e-12 for a, b in zip(e["min"], e["ece"]))
    assert all(a <= b + 1e-12 for a, b in zip(e["min"], e["ref"]))
    assert len(eceplt(s, d, draw=False)["lo"]) == 101


def test_eceplt_takes_the_limits_at_extreme_odds():
    e = eceplt([2, 3], 0.5, lo=[-400, 0, 400], draw=False)
    assert e["min"] == [0.0, 0.0, 0.0]
    want = (1e-100 * 400 * math.log(10) + 1e-100) / math.log(2)
    assert eceplt(1e-300, 1, lo=-100, draw=False)["ece"][0] == pytest.approx(want, rel=1e-12)
    want = 0.5 * (310 * math.log2(10) + 1)
    assert eceplt(1e-310, 1, lo=0, draw=False)["ece"][0] == pytest.approx(want, rel=1e-12)


def test_evidence_weights_never_underflow():
    p = [1e-3] * 500
    g = [[i, i] for i in range(500)]
    r = kinidx(g, g, [[v, v] for v in p], rel=[0.9999, 0.0001, 0])
    assert r["log10idx"] == pytest.approx(500 * math.log10(0.9999 + 0.0001 / 1e-3), rel=1e-10)
    assert lrpost(1e-300, 1 - 1e-16)[0] == pytest.approx(1 / (1 + 1e-16 / (1 - 1e-16) * 1e300),
                                                          rel=1e-10)


def test_frequencies_outside_their_range_are_refused():
    with pytest.raises(ValueError, match="at most 1"):
        mixlr({"a": 0.5, "b": 0.5 + 1e-9}, "a", 1, "a", 1)
    with pytest.raises(ValueError, match="`q`"):
        gtmp(0.5, 0.5 + 1e-9)
    with pytest.raises(ValueError, match="`freq`"):
        kinidx([1, 2], [1, 2], [0, 0.2])
    with pytest.raises(ValueError, match="`freq`"):
        patidx([1, 2], [1, 2], [0.1, 1.01], mother=[1, 1])
    with pytest.raises(TypeError, match="`freq`"):
        kinidx([1, 2], [1, 2], ["a", "b"])


def test_a_group_never_seen_in_daylight_is_no_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        vodtst({"a": 50, "b": 0}, {"a": 40, "b": 3})


def test_orrisk_refuses_a_reciprocal_that_overflows():
    with pytest.raises(ValueError, match="reciprocal"):
        orrisk(1e-320, 0.5, 0.5)


def test_number_checks():
    assert incap([1, 2], 0.5, 1)["lam"] == [1.0, 2.0]
    assert incap(1, 0.5, 0)["s"] == [0.0]
    assert hzdsel(0.5, 1, 0, [1, 1], [1, 1])["share0"] == 0.5
    with pytest.raises(ValueError, match="length 3, not 2"):
        rtmhs([1, 2, 3], [1, 2], 0)
    with pytest.raises(ValueError, match="at least 0"):
        losbnd(-1, 1)
    with pytest.raises(ValueError, match="greater than 0"):
        cllr(0, 1)
    with pytest.raises(ValueError, match="strictly between 0 and 1"):
        orrisk(2, 0, 0.5)


def test_the_reference_by_name_or_position():
    assert expbnd({"a": 1, "b": 2}, [1, 1], 1, ref="b")["ratio"] == [0.5, 1.0]
    assert expbnd({"a": 1, "b": 2}, [1, 1], 1, ref=2)["ratio"] == [0.5, 1.0]
    assert vodtst([4, 5], [1, 1])["group"] == ["1", "2"]
    assert vodtst({"a": 4, "b": 5}, [1, 1])["group"] == ["a", "b"]
    with pytest.raises(ValueError, match="position from 1 to 2"):
        vodtst([1, 1], [1, 1], ref=1.5)
    with pytest.raises(TypeError, match="`ref`"):
        vodtst([1, 1], [1, 1], ref=True)
    with pytest.raises(TypeError, match="`ref`"):
        vodtst([1, 1], [1, 1], ref=None)
    with pytest.raises(ValueError, match="names no group: z"):
        vodtst({"a": 1, "b": 1}, [1, 1], ref="z")


def test_a_genotype_of_two_values_is_one_row():
    r = kinidx([1, 2], [1, 2], [0.1, 0.2])
    assert kinidx(["1", "2"], [1, 2], [0.1, 0.2]) == r
    assert kinidx([1.0, 2.0], ("1", "2"), [[0.1, 0.2]]) == r
    assert kinidx([0.1 + 0.2, 0.3], [0.3, 0.3], [0.1, 0.1]) == kinidx([0.3, 0.3], [0.3, 0.3],
                                                                      [0.1, 0.1])
    with pytest.raises(TypeError, match="`g1`"):
        kinidx([True, False], [1, 2], [0.1, 0.2])
    with pytest.raises(TypeError, match="`g1`"):
        kinidx([1, object()], [1, 2], [0.1, 0.2])
    with pytest.raises(ValueError, match="`g1`"):
        kinidx("ab", [1, 2], [0.1, 0.2])


def test_many_alleles_are_coded_alike_by_value_and_by_label():
    g = [[i, 3000 + i] for i in range(1, 3001)]
    h = [[3000 + i, i] for i in range(1, 3001)]
    f = [[0.001, 0.001]] * 3000

    def t(x):
        return [[str(a) for a in r] for r in x]

    assert kinidx(g, h, f) == kinidx(t(g), t(h), f)
    assert patidx(g, h, f, g) == patidx(t(g), t(h), f, t(g))


def test_t0_plus_s_past_the_largest_integer():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match=r"length at least t0 \+ s \+ 1"):
            incapr([2, 1], 2**31 - 1, 1)


def test_z_is_two_labels_and_treat_is_never_chosen_by_sorting():
    assert sntbnd([0, 1], ["b", "a"], "b")["term"] == ["a", "b", "effect"]
    assert sntbnd([0, 1], ["a", "B"], "a")["term"] == ["B", "a", "effect"]
    assert sntbnd([0, 1], [2, 10], 10)["term"] == ["2", "10", "effect"]
    assert sntbnd([0, 1], [2.0, 10.0], 10.0) == sntbnd([0, 1], [2, 10], "10")
    assert sntbnd([0, 1], ("a", "b"), "b") == sntbnd([0, 1], ["a", "b"], "b")
    assert sntbnd([0, 1], [True, False])["term"] == ["FALSE", "TRUE", "effect"]
    assert sntbnd([1, 0, 1, 0], [1, 0, 1, 0])["term"] == ["0", "1", "effect"]
    assert sntbnd([1, 0, 1, 0], [1, 0, 1, 0]) == sntbnd([1, 0, 1, 0], [1, 0, 1, 0], 1)
    with pytest.raises(ValueError, match="`treat` must name the treatment, a or b"):
        sntbnd([0, 1], ["a", "b"])
    with pytest.raises(ValueError, match="`treat` must name the treatment, 1 or 2"):
        sntbnd([0, 1], [1, 2])
    with pytest.raises(ValueError, match="`treat` must name the treatment, 0 or 1"):
        sntbnd([0, 1], ["0", "1"])
    with pytest.raises(ValueError, match="exactly two"):
        sntbnd([0, 1, 1], ["a", "b", "c"])
    with pytest.raises(TypeError, match="labels"):
        sntbnd([0, 1], [object(), object()])
    # a mapping is labels named by unit
    assert sntbnd({"u": 0, "v": 1}, {"v": 2, "u": 1}, 2) == sntbnd([0, 1], [1, 2], 2)
    with pytest.raises(ValueError, match="`treat`"):
        sntbnd([0, 1], ["a", "b"], treat=["a"])
    with pytest.raises(ValueError, match="`treat`"):
        sntbnd([0, 1], ["a", "b"], treat=math.nan)
    with pytest.raises(ValueError, match="`treat`"):
        sntbnd([0, 1], ["a", "b"], treat=object())
    with pytest.raises(ValueError, match="no missing"):
        sntbnd([0, 1], ["a", math.nan])
    assert sntbnd([0, 1], ["", "a"], treat="")["term"] == ["a", "", "effect"]
