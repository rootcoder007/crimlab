# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""NumPy, pandas and matplotlib, when they are installed (the package
imports none of them): NumPy arrays as inputs through their array
interface, zero copy where contiguous float64; pandas Series and
DataFrames as mappings through items(), keeping their labels; results turned into theirs by their own
constructors; the plot data drawn with matplotlib. The only test file that
imports them. Not collected by
default (its name does not start with test_): run it as
`python -m pytest opt.py` where they are installed; skipped without."""

import pytest

from crimlab import cllr, dyadrt, gtmp, losbnd, sellab, sntbnd, vodtst

DONE = [30.0, 45.0, 60.0, 90.0, 120.0]

np = pytest.importorskip("numpy")


def test_numpy_arrays_of_every_numeric_type():
    want = losbnd(DONE, 100)
    for t in ("float64", "float32", "int64", "int32", "int16", "int8", "uint8", "uint16",
              "uint32", "uint64", "float16"):
        assert losbnd(np.array(DONE).astype(t), 100) == want, t
    assert losbnd(np.array(DONE).astype(">f8"), 100) == want
    assert losbnd(np.array(DONE, dtype=object), 100) == want
    assert losbnd(np.float64(30), 100) == losbnd(30, 100)
    assert losbnd(np.array(30.0), np.int64(100)) == losbnd(30, 100)


def test_numpy_arrays_not_contiguous_and_of_more_dimensions():
    x = np.arange(20, dtype=float) + 1
    assert losbnd(x[::2], 5) == losbnd(list(x[::2]), 5)
    assert losbnd(x[::-1], 5) == losbnd(list(x[::-1]), 5)
    assert losbnd(x[::3].astype("int32"), 5) == losbnd(list(x[::3]), 5)
    assert losbnd(x.reshape(20, 1), 5) == losbnd(x, 5)
    assert losbnd(x.reshape(1, 20), 5) == losbnd(x, 5)
    assert losbnd(x.reshape(1, 20, 1), 5) == losbnd(x, 5)
    assert losbnd(np.asfortranarray(x.reshape(20, 1)), 5) == losbnd(x, 5)
    m = "must be a vector, not a matrix with more than one row and column"
    with pytest.raises(ValueError, match=m):
        losbnd(x.reshape(4, 5), 5)
    with pytest.raises(ValueError, match=m):
        vodtst(np.array([[6.0, 4.0], [5.0, 5.0]]), [5, 5, 4, 6])
    with pytest.raises(ValueError, match=m):
        losbnd(np.array([[1, None], [3, 4]], dtype=object), 5)
    with pytest.raises(ValueError, match=m):
        sellab([1, 0, 1, 0], np.array([[True, False], [True, True]]), [True] * 4)
    with pytest.raises(ValueError, match=m):
        sellab([1, 0, 1, 0], np.array([[True, False], [True, True]], dtype=object), [True] * 4)
    with pytest.raises(ValueError, match=m):
        sellab(np.ones((2, 2)), [True] * 4, [True] * 4)
    with pytest.raises(ValueError, match=m):
        gtmp([0.1, 0.2, 0.1, 0.2], np.ones((2, 2)) / 20)
    with pytest.raises(ValueError, match=m):
        gtmp([0.1, 0.2, 0.1, 0.2], hom=np.ones((2, 2), dtype=bool))


def test_numpy_arrays_that_are_not_numbers():
    for bad in (np.array([True, False]), np.array(["1", "2"]), np.array([1 + 1j, 2]),
                np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[D]")):
        with pytest.raises(TypeError, match="`done` must be finite numbers"):
            losbnd(bad, 1)
    with pytest.raises(ValueError, match="`done` must be finite numbers"):
        losbnd(np.array([1.0, np.nan]), 1)
    assert sellab(np.array([1.0, 0.0]), np.array([True, True]), np.array([True, False])) == \
        sellab([1, 0], [True, True], [True, False])
    with pytest.raises(TypeError, match="`rel` must be TRUE or FALSE"):
        sellab([1, 0], np.array([1, 0]), [True, True])
    assert sellab([1, 0], np.array([True, True], dtype=object), [np.bool_(True), True]) == \
        sellab([1, 0], [True, True], [True, True])
    with pytest.raises(TypeError, match="`done` must be finite numbers"):
        losbnd([np.bool_(True)], 1)


def test_numpy_results_and_inputs_agree_with_lists():
    rng = np.random.default_rng(1)
    s = rng.lognormal(1, 1, 1000)
    d = rng.lognormal(-1, 1, 1000)
    assert cllr(s, d) == cllr(s.tolist(), d.tolist())


def test_pandas_objects_keep_their_labels():
    pd = pytest.importorskip("pandas")
    # a Series is a mapping through items(): two in different label order
    # are matched by name, never by position
    day = pd.Series({"w": 600.0, "b": 400.0})
    night = pd.Series({"b": 300.0, "w": 700.0})
    assert vodtst(day, night) == vodtst({"w": 600, "b": 400}, {"w": 700, "b": 300})
    with pytest.raises(ValueError, match="`day` names a group twice: a"):
        vodtst(pd.Series([1.0, 2.0], index=["a", "a"]), [1, 1])
    # a DataFrame is a mapping of columns, each by row label
    t = pd.DataFrame({"night": [700, 300], "day": [600, 400]}, index=["w", "b"])
    assert vodtst(t) == vodtst({"w": 600, "b": 400}, {"w": 700, "b": 300})
    c = pd.DataFrame([[640, 160], [100, 40]], index=["a", "b"], columns=["a", "b"])
    pop = {"a": 800, "b": 200}
    assert dyadrt(c.loc[["b", "a"], ["b", "a"]], pop) == dyadrt([[640, 160], [100, 40]], pop)
    with pytest.raises(TypeError, match="`done` must be finite numbers"):
        losbnd([30, pd.NA], 100)
    with pytest.raises(TypeError, match="`night` must be finite numbers"):
        vodtst(day, pd.Series({"w": 700, "b": pd.NA}, dtype="object"))
    # a categorical is read by its labels
    y = [0, 1, 1, 0]
    assert sntbnd(y, pd.Categorical(["x", "y", "x", "y"], categories=["y", "x"]), "y") == \
        sntbnd(y, ["x", "y", "x", "y"], "y")


def test_a_result_is_a_data_frame_by_pandas_own_constructor():
    pd = pytest.importorskip("pandas")
    r = vodtst({"w": 600, "b": 400}, {"w": 700, "b": 300})
    f = pd.DataFrame(r)
    assert list(f.columns) == list(r) and f.shape == r.shape
    assert pd.DataFrame(r.filter(ref=False))["group"].tolist() == ["b"]


def test_a_column_is_an_array_by_numpy_own_constructor():
    np = pytest.importorskip("numpy")
    r = vodtst({"w": 600, "b": 400}, {"w": 700, "b": 300})
    a = np.asarray(r["ratio"])
    assert a.dtype == np.float64 and a.tolist() == r["ratio"]


def test_the_plot_data_draws_with_matplotlib():
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    from crimlab import tipplt

    t = tipplt([50, 200, 8, 0.5], [0.01, 0.2, 1.5, 0.001])
    fig, ax = plt.subplots()
    for k in ("same", "diff"):
        ax.step(t[k]["log10lr"], t[k]["share"], where="post", label=k)
    assert len(ax.lines) == 2
    plt.close(fig)
