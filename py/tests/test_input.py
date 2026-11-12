# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""What the binding accepts as numbers: each kind of input has a stated
behaviour (see README.md, "Inputs")."""

import array
import math
from decimal import Decimal
from fractions import Fraction

import pytest

from crimlab import cllr, gtmp, hzdsel, losbnd, lrpost, outtst, sellab, vodtst

DONE = [30.0, 45.0, 60.0, 90.0, 120.0]
WANT = losbnd(DONE, [100, 150, 200])


@pytest.mark.parametrize("x", [
    DONE, tuple(DONE), array.array("d", DONE), array.array("i", [30, 45, 60, 90, 120]),
    memoryview(array.array("d", DONE)), [30, 45, 60, 90, 120], range(30, 31),
    [Decimal(v) for v in (30, 45, 60, 90, 120)], [Fraction(v) for v in DONE],
])
def test_sequences_and_buffers(x):
    if isinstance(x, range):
        assert losbnd(x, [100, 150, 200]) == losbnd([30], [100, 150, 200])
        return
    assert losbnd(x, [100, 150, 200]) == WANT


def test_one_number_is_a_vector_of_one():
    assert losbnd(DONE, 100) == losbnd(DONE, [100])
    assert lrpost(2, 0.5) == lrpost([2], [0.5])


@pytest.mark.parametrize("x, err", [
    ("abc", TypeError), (b"ab", TypeError), (None, TypeError), ({"a": 1}, TypeError),
    ({1, 2}, TypeError), ((v for v in DONE), TypeError), ([[30, 45]], TypeError),
    ([True, False], TypeError), (True, TypeError), (["1"], TypeError), ([1, object()], TypeError),
    ([], ValueError), ([1, None], ValueError), ([1, math.nan], ValueError),
    ([1, math.inf], ValueError), ([10**400], ValueError), ([-(10**400)], ValueError),
])
def test_what_is_not_numbers(x, err):
    with pytest.raises(err, match=r"^losbnd: `done` must be finite numbers with no missing values\.$"):
        losbnd(x, 1)


def test_scalar_arguments():
    want = hzdsel(0.5, 0.5, 0.1, [0.5, 0.8], [0.9, 0.9])
    assert hzdsel([0.5], (0.5,), Decimal("0.1"), [0.5, 0.8], [0.9, 0.9]) == want
    with pytest.raises(TypeError, match="`s` must be one number"):
        hzdsel("0.5", 0.5, 0.1, [0.5, 0.8], [0.9, 0.9])
    with pytest.raises(ValueError, match="`s` must be one number"):
        hzdsel([0.5, 0.5], 0.5, 0.1, [0.5, 0.8], [0.9, 0.9])
    with pytest.raises(ValueError, match="`s` must be one number"):
        hzdsel(math.nan, 0.5, 0.1, [0.5, 0.8], [0.9, 0.9])
    with pytest.raises(TypeError, match="`s` must be one number"):
        hzdsel(None, 0.5, 0.1, [0.5, 0.8], [0.9, 0.9])


def test_logical_arguments():
    want = sellab([1, 0, None], [True, True, False], [True, False, True])
    assert sellab([1, 0, None], (True, True, False), [True, False, True]) == want
    with pytest.raises(TypeError, match="`rel` must be TRUE or FALSE"):
        sellab([1, 0], [1, 0], [True, True])
    with pytest.raises(ValueError, match="`rel` must be TRUE or FALSE"):
        sellab([1, 0], [True, None], [True, True])
    with pytest.raises(TypeError, match="`rel` must be TRUE or FALSE"):
        sellab([1, 0], "TT", [True, True])
    with pytest.raises(TypeError, match="`hom`"):
        gtmp(0.2, hom=1)


def test_very_long_inputs():
    n = 1_000_000
    s = [1.0 + (i % 7) for i in range(n)]
    d = [1.0 / (1 + (i % 5)) for i in range(n)]
    assert cllr(s, d) == cllr(tuple(s), array.array("d", d))
    assert outtst([0.5] * n, 0.3) == {"stop": 1.0, "hit": 0.5}


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


pd = pytest.importorskip("pandas")


def test_pandas_series_and_missing_values():
    assert losbnd(pd.Series(DONE), 100) == losbnd(DONE, 100)
    assert losbnd(pd.Series([30, 45, 60, 90, 120], dtype="Int64"), 100) == losbnd(DONE, 100)
    with pytest.raises(ValueError, match="`done` must be finite numbers"):
        losbnd(pd.Series([30, None], dtype="Float64"), 100)
    with pytest.raises(ValueError, match="`done` must be finite numbers"):
        losbnd([30, pd.NA], 100)
    with pytest.raises(TypeError, match="`done` must be finite numbers"):
        losbnd(pd.Series(["a", "b"]), 100)
    with pytest.raises(TypeError, match="`done` must be finite numbers"):
        losbnd(pd.Series([True, False]), 100)
    assert sellab(pd.Series([1, 0]), pd.Series([True, True]), pd.Series([True, False])) == \
        sellab([1, 0], [True, True], [True, False])
