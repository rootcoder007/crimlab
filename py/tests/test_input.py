# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""What the binding accepts as numbers: each kind of input has a stated
behaviour (see README.md, "Inputs")."""

import array
import math
import sys
from decimal import Decimal
from fractions import Fraction

import pytest

from crimlab import cllr, gtmp, hzdsel, losbnd, lrpost, outtst, sellab

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


class Arr:
    """An array over memory the test owns: what an array library offers
    through __array_interface__, built from array.array."""

    def __init__(self, code, values, typestr, shape=None, strides=None):
        self.buf = array.array(code, values)
        addr, n = self.buf.buffer_info()
        self.__array_interface__ = {"typestr": typestr, "shape": shape or (n,),
                                    "strides": strides, "data": (addr, True), "version": 3}


E = "<" if sys.byteorder == "little" else ">"


def test_array_interfaces_of_every_numeric_type():
    want = losbnd(DONE, 100)
    for code, ts in (("d", "f8"), ("f", "f4"), ("b", "i1"), ("h", "i2"), ("i", "i4"),
                     ("q", "i8"), ("B", "u1"), ("H", "u2"), ("I", "u4"), ("Q", "u8")):
        assert losbnd(Arr(code, [int(v) for v in DONE] if code not in "df" else DONE,
                          E + ts), 100) == want, ts


def test_array_interfaces_strided_and_of_more_dimensions():
    x = [float(v) for v in range(1, 21)]
    assert losbnd(Arr("d", x, E + "f8", (10,), (16,)), 5) == losbnd(x[::2], 5)
    assert losbnd(Arr("i", list(range(1, 21)), E + "i4", (10,), (8,)), 5) == \
        losbnd(x[::2], 5)
    assert losbnd(Arr("d", x, E + "f8", (20, 1)), 5) == losbnd(x, 5)
    assert losbnd(Arr("d", x, E + "f8", (1, 20, 1)), 5) == losbnd(x, 5)
    m = "must be a vector, not a matrix with more than one row and column"
    with pytest.raises(ValueError, match=m):
        losbnd(Arr("d", x, E + "f8", (4, 5)), 5)
    with pytest.raises(ValueError, match=m):
        sellab([1, 0, 1, 0], Arr("b", [1, 0, 1, 1], "|b1", (2, 2)), [True] * 4)
    with pytest.raises(ValueError, match=m):
        gtmp([0.1, 0.2, 0.1, 0.2], Arr("d", [0.05] * 4, E + "f8", (2, 2)))
    # a byte order not the machine's is not read from memory: an object
    # that is no sequence as well is then not numbers
    with pytest.raises(TypeError, match="`done` must be finite numbers"):
        losbnd(Arr("d", x, ("<" if E == ">" else ">") + "f8"), 5)


def test_array_interfaces_that_are_not_numbers():
    for ts in ("|b1", E + "c16", E + "M8[D]", "<U1"):
        with pytest.raises(TypeError, match="`done` must be finite numbers"):
            losbnd(Arr("d", [1.0, 2.0], ts), 1)
    assert sellab([1, 0], Arr("b", [1, 1], "|b1"), Arr("b", [1, 0], "|b1")) == \
        sellab([1, 0], [True, True], [True, False])
    with pytest.raises(TypeError, match="`rel` must be TRUE or FALSE"):
        sellab([1, 0], Arr("i", [1, 0], E + "i4"), [True, True])
