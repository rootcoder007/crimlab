# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""The plots as SVG of our own (src/svg.cpp, clab_svgln)."""

import math
import re
import xml.etree.ElementTree as ET

import pytest

import crimlab
from crimlab import _clab

NS = "{http://www.w3.org/2000/svg}"
S = [50, 200, 8, 0.5]
D = [0.01, 0.2, 1.5, 0.001]


def curves(svg):
    root = ET.fromstring(svg.encode())
    return root, [p for p in root.iter(NS + "path") if "curve" in p.get("class", "")]


def first(d):
    m = re.match(r"M([-0-9.]+) ([-0-9.]+)", d)
    return float(m.group(1)), float(m.group(2))


def last(d):
    nums = re.findall(r"([HVL ])([-0-9.]+)", d)
    x = y = None
    for k, v in nums:
        if k == "H":
            x = float(v)
        elif k == "V":
            y = float(v)
    return x, y


def px(v, lo, hi, a, b):
    p = 0.04 * (hi - lo)
    lo, hi = lo - p, hi + p
    return a + (v - lo) / (hi - lo) * (b - a)


def test_tippett_svg():
    t = crimlab.tipplt(S, D)
    s = t.svg()
    root, cs = curves(s)
    assert root.tag == NS + "svg" and root.get("viewBox") == "0 0 640 400"
    assert len(cs) == 2
    assert root.find(NS + "title").text == "Tippett plot"
    texts = [e.text for e in root.iter(NS + "text")]
    assert "same source" in texts and "different source" in texts
    assert "log10 likelihood ratio" in texts
    # the first point of the same-source curve and the last of the other
    xs = t["same"]["log10lr"] + t["diff"]["log10lr"] + [0]
    lo, hi = min(xs), max(xs)
    x, y = first(cs[0].get("d"))
    assert x == pytest.approx(px(t["same"]["log10lr"][0], lo, hi, 64, 624), abs=0.006)
    assert y == pytest.approx(px(1.0, 0, 1, 348, 40), abs=0.006)
    x, y = last(cs[1].get("d"))
    assert x == pytest.approx(px(t["diff"]["log10lr"][-1], lo, hi, 64, 624), abs=0.006)
    assert y == pytest.approx(px(0.25, 0, 1, 348, 40), abs=0.006)


def test_ece_svg():
    e = crimlab.eceplt(S, D)
    root, cs = curves(e.svg(title="ECE & <more>"))
    assert len(cs) == 3
    assert root.find(NS + "title").text == "ECE & <more>"
    x, y = first(cs[2].get("d"))
    hi = max(e["ece"] + e["min"] + e["ref"])
    assert x == pytest.approx(64.0 + 0.04 / 1.08 * 560, abs=0.006)
    assert y == pytest.approx(px(e["ref"][0], 0, hi, 348, 40), abs=0.006)


def test_svg_bytes_repeat_and_text_is_small():
    a = crimlab.tipplt(S, D).svg()
    assert a == crimlab.tipplt(S, D).svg() == crimlab.tipplt(list(S), tuple(D)).svg()
    assert "font-family:system-ui" in a and "prefers-color-scheme:dark" in a
    assert len(a.encode()) < 8000


def test_special_values():
    # a single ratio on each side: a valid picture
    cs = curves(crimlab.tipplt([3], [0.5]).svg())[1]
    assert len(cs) == 2
    # one lo: a single x is widened
    cs = curves(crimlab.eceplt([2, 3], [0.5], lo=[0]).svg())[1]
    assert len(cs) == 3
    # extreme odds: calibrated cost 0, still finite
    curves(crimlab.eceplt([2, 3], [0.5], lo=[-400, 0, 400]).svg())
    # no ratios on a side is refused by tipplt itself, as in R
    with pytest.raises(ValueError, match="`same`"):
        crimlab.tipplt([], [1])
    # infinities never reach a picture: a clear error
    with pytest.raises(ValueError, match="outside the domain"):
        _clab.svg([([0.0, 1.0], [0.0, math.inf], False, 0, "x")], "t", "x", "y", False, 0.0, 0)
    with pytest.raises(ValueError, match="out of range"):
        _clab.svg([([0.0], [0.0], False, 0, "x")] * 4, "t", "x", "y", False, 0.0, 0)
    with pytest.raises(ValueError, match="out of range"):
        _clab.svg([([0.0], [0.0], False, 0, "x")], "t", "x", "y", False, 0.0, 4)
    with pytest.raises(TypeError, match="a curve must be two lists of numbers of one length"):
        _clab.svg([([0.0], [0.0, 1.0], False, 0, "x")], "t", "x", "y", False, 0.0, 0)
    with pytest.raises(TypeError, match="a curve must be"):
        _clab.svg([([0.0], [0.0])], "t", "x", "y", False, 0.0, 0)
    with pytest.raises(ValueError, match="dash must be 0, 1 or 2"):
        _clab.svg([([0.0], [0.0], False, 3, "x")], "t", "x", "y", False, 0.0, 0)
