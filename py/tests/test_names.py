# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Names are the identity of a group: two named inputs are matched by
name, the output keeps the order of the first, and the reference group is
never chosen by sorting."""

import pytest

from crimlab import (
    clrsmp,
    dspbm,
    dyadrt,
    expbnd,
    gtmp,
    incap,
    kinidx,
    mixlr,
    patidx,
    sntbnd,
    vodtst,
)


def pd_series(v, idx):
    pd = pytest.importorskip("pandas")
    return pd.Series(v, index=idx)


def rev(d):
    return dict(reversed(list(d.items())))


DAY = {"w": 600, "b": 400, "h": 150}
NIGHT = {"w": 700, "b": 300, "h": 140}
SHIFT = {"w": 1.0, "b": 1.2, "h": 0.9}


def test_named_inputs_are_matched_by_name():
    want = vodtst(list(DAY.values()), list(NIGHT.values()), shift=list(SHIFT.values()))
    got = vodtst(DAY, rev(NIGHT), shift=rev(SHIFT))
    assert got["group"] == ["w", "b", "h"]
    assert {k: v for k, v in got.items() if k != "group"} == \
        {k: v for k, v in want.items() if k != "group"}


@pytest.mark.parametrize("call", [
    lambda a, b, c: dspbm(a, b, c),
    lambda a, b, c: clrsmp(a, b, a, c),
    lambda a, b, c: expbnd(a, b, 1.5),
    lambda a, b, c: incap(a, {k: 0.1 for k in b}, 1, share={k: 1 / 3 for k in c}),
])
def test_every_order_of_the_later_names_gives_one_answer(call):
    a = {"x": 9.0, "y": 3.0, "z": 5.0}
    b = {"x": 3.0, "y": 2.0, "z": 4.0}
    c = {"x": 2.0, "y": 1.0, "z": 7.0}
    ref = call(a, b, c)
    for perm in (["y", "z", "x"], ["z", "x", "y"], ["z", "y", "x"]):
        assert call(a, {k: b[k] for k in perm}, {k: c[k] for k in perm}) == ref


def test_the_first_input_orders_the_output_and_picks_the_reference():
    r = vodtst({"b": 4, "a": 6}, {"a": 5, "b": 5})
    assert r["group"] == ["b", "a"]
    assert r["ratio"][0] == 1.0
    r = expbnd({"z": 1, "a": 2}, {"a": 1, "z": 1}, 1)
    assert r["group"] == ["z", "a"] and r["ratio"] == [1.0, 2.0]


def test_different_names_are_an_error_listing_them():
    with pytest.raises(ValueError) as e:
        vodtst({"a": 1, "b": 2}, {"a": 1, "c": 2})
    assert str(e.value) == "vodtst: `night` must name the same groups as `day`; missing: b; extra: c."
    with pytest.raises(ValueError, match="`q` must name the same groups as `p`"):
        gtmp({"D3": 0.1, "TH01": 0.2}, {"vWA": 0.2, "TH01": 0.3})
    with pytest.raises(ValueError, match="`day` names a group twice: a"):
        vodtst(pd_series([1, 2], ["a", "a"]), [1, 1])
    # a wrong length is not matched by name: the length check refuses it
    with pytest.raises(ValueError, match="`contact` must have length 3, not 2"):
        dspbm({"a": 1, "b": 1, "c": 1}, {"a": 1, "b": 1}, [1, 1, 1])


def test_an_unnamed_input_is_taken_by_position():
    assert vodtst({"a": 6, "b": 4}, [5, 5]) == vodtst({"a": 6, "b": 4}, {"a": 5, "b": 5})
    # the groups are named by the first argument that names them
    assert vodtst([6, 4], {"y": 5, "x": 5})["group"] == ["y", "x"]
    assert vodtst([6, 4], [5, 5])["group"] == ["1", "2"]
    # one table of day and night
    assert vodtst({"night": {"a": 5, "b": 5}, "day": {"a": 6, "b": 4}}) == \
        vodtst({"a": 6, "b": 4}, {"a": 5, "b": 5})
    with pytest.raises(ValueError, match="`night` is missing"):
        vodtst({"a": 6, "b": 4})


def test_per_locus_arguments_are_matched_by_locus():
    a = gtmp({"D3": 0.1, "vWA": 0.2}, {"vWA": None, "D3": 0.3}, hom={"vWA": True, "D3": False})
    b = gtmp({"D3": 0.1, "vWA": 0.2}, [0.3, None], hom=[False, True])
    assert a == b


def test_genotype_tables_are_matched_by_locus():
    g1 = {"D3": [1, 2], "vWA": [3, 3]}
    g2 = {"vWA": [3, 4], "D3": [2, 2]}
    f = {"vWA": [0.3, 0.3], "D3": [0.1, 0.2]}
    a = kinidx(g1, g2, f, "sibling")
    b = kinidx([[1, 2], [3, 3]], [[2, 2], [3, 4]], [[0.1, 0.2], [0.3, 0.3]], "sibling")
    assert a["locus"]["locus"] == ["D3", "vWA"]
    assert a["locus"]["idx"] == b["locus"]["idx"]
    m = {"vWA": [3, 6], "D3": [1, 1]}
    assert patidx(g1, g2, f, m)["locus"]["idx"] == \
        patidx([[1, 2], [3, 3]], [[2, 2], [3, 4]], [[0.1, 0.2], [0.3, 0.3]],
               [[1, 1], [3, 6]])["locus"]["idx"]
    with pytest.raises(ValueError, match="`g2` must name the same groups as `g1`"):
        kinidx(g1, {"D3": [1, 1], "TH01": [2, 2]}, f)


def test_mixture_alleles_are_names():
    p = {"a": 0.1, "b": 0.2}
    assert mixlr(p, "a", 1, ("b",), 1) == mixlr(rev(p), ["a"], 1, ["b"], 1)
    assert mixlr(p, {"a"}, 1, {"b"}, 1) == mixlr(p, "a", 1, "b", 1)
    with pytest.raises(ValueError, match="uniquely named"):
        mixlr([0.1, 0.2], "a", 1, "a", 1)
    with pytest.raises(ValueError, match="`np` must be names"):
        mixlr(p, [["a"]], 1, "a", 1)
    with pytest.raises(ValueError, match="`np` must be names"):
        mixlr(p, 3, 1, "a", 1)
    with pytest.raises(ValueError, match="at most 20"):
        mixlr({str(i): 0.01 for i in range(21)}, [str(i) for i in range(21)], 1, "0", 1)


pd = pytest.importorskip("pandas")


def test_series_are_matched_by_their_index():
    day = pd.Series(DAY)
    night = pd.Series(NIGHT)[["h", "w", "b"]]
    assert vodtst(day, night) == vodtst(DAY, NIGHT)
    # a RangeIndex is positions, not names
    assert vodtst(pd.Series([6, 4]), pd.Series([5, 5]))["group"] == ["1", "2"]
    with pytest.raises(ValueError, match="`day` names a group twice: a"):
        vodtst(pd.Series([1, 2], index=["a", "a"]), pd.Series([1, 2], index=["a", "b"]))


def test_a_categorical_is_read_by_its_labels():
    y = [0, 1, 1, 0]
    z = ["x", "y", "x", "y"]
    want = sntbnd(y, z, "y")
    assert sntbnd(y, pd.Categorical(z, categories=["y", "x"]), "y") == want
    assert sntbnd(y, pd.Series(z, dtype="category"), "y") == want
    with pytest.raises(ValueError, match="`treat` must name"):
        sntbnd(y, pd.Categorical([0, 1, 0, 1]))
    with pytest.raises(TypeError, match="`day` must be finite numbers"):
        vodtst(pd.Series([1, 2], dtype="category"), [1, 1])


def test_a_labelled_count_table_follows_pop():
    pop = {"a": 800, "b": 200}
    c = pd.DataFrame([[640, 160], [100, 40]], index=["a", "b"], columns=["a", "b"])
    want = dyadrt(c.to_numpy(), pop)
    assert dyadrt(c.loc[["b", "a"], ["b", "a"]], pop) == want
    with pytest.raises(ValueError, match="the row names of `count` must name the same groups"):
        dyadrt(c.rename(index={"b": "c"}), pop)


def test_genotypes_and_frequencies_as_data_frames():
    g = pd.DataFrame({"a1": [1, 3], "a2": [2, 3]}, index=["D3", "vWA"])
    f = pd.DataFrame({"a1": [0.1, 0.3], "a2": [0.2, 0.3]}, index=["D3", "vWA"])
    r = kinidx(g, g.iloc[::-1], f.iloc[::-1])
    assert r == kinidx({"D3": [1, 2], "vWA": [3, 3]}, {"D3": [1, 2], "vWA": [3, 3]},
                       {"D3": [0.1, 0.2], "vWA": [0.3, 0.3]})
    cat = g.astype("category")
    assert kinidx(cat, g, f) == kinidx(g, g, f)
