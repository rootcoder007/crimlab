# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Memory and references: every call, valid or not, releases what it took,
and a failed allocation in the compiled module is a MemoryError, never a
crash, a wrong result or a leak. The failing allocations need CPython's _testcapi and
CLAB_NOMEM=1 (they take minutes); the reference counts a debug build
(sys.gettotalrefcount); each test is skipped where its tool is missing."""

import gc
import os
import subprocess
import sys
from pathlib import Path

import pytest
from ref import CASES, args

import crimlab
from crimlab import _clab

HERE = Path(__file__).parent


def run(k, fn):
    for _ in range(k):
        for c in CASES:
            if c["fn"] == fn:
                try:
                    getattr(crimlab, fn)(**args(c))
                except (ValueError, TypeError):
                    pass


@pytest.mark.skipif(not hasattr(sys, "gettotalrefcount"), reason="needs a debug build")
def test_reference_counts_do_not_grow():
    # A leak adds the same count in every window; the interpreter's own
    # caches add some at odd times while they fill (slowly on 3.15). So per
    # function, after every function has run: the least growth over three
    # windows of 30 rounds of its cases, against the number of calls in a
    # window, the least a leak of one reference per call would add.
    fns = sorted({c["fn"] for c in CASES})
    for fn in fns:
        run(5, fn)
    for fn in fns:
        grow = []
        for _ in range(3):
            gc.collect()
            a = sys.gettotalrefcount()
            run(30, fn)
            gc.collect()
            grow.append(sys.gettotalrefcount() - a)
        calls = 30 * sum(c["fn"] == fn for c in CASES)
        assert min(grow) < calls, (fn, grow, calls)


NOMEM = """
import json, sys, _testcapi, crimlab
from crimlab import _clab
from ref import CASES, args
fail, unhook = _testcapi.set_nomemory, _testcapi.remove_mem_hooks
# The calls each wrapper makes into the compiled module are recorded, then
# replayed with one allocation failing at a time: only the extension runs
# while an allocation can fail. (Failing allocations inside Python bytecode
# itself aborts CPython 3.15's interpreter, with or without crimlab.)
real = {n: getattr(_clab, n) for n in dir(_clab) if not n.startswith("__")}
seen = []
for name, fn in real.items():
    def spy(*a, _n=name, _f=fn):
        seen.append((_n, a))
        return _f(*a)
    setattr(_clab, name, spy)
n = 0
# the small cases reach every allocation site; a large one would fail each
# of its thousands of allocations in turn, at a cost quadratic in its size
for i, c in enumerate(CASES):
    if len(json.dumps(c["args"])) >= 800:
        continue
    seen.clear()
    try:
        getattr(crimlab, c["fn"])(**args(c))
    except (ValueError, TypeError):
        pass
    for name, a in list(seen):
        f = real[name]
        try:
            want = repr(f(*a))
        except (ValueError, TypeError) as e:
            want = repr(e)
        # the case and the call go to stderr first, so that a crash still
        # names them in the log
        sys.stderr.write(f"case {i} ({c['fn']}): {name}\\n")
        sys.stderr.flush()
        # one allocation fails at a time: the start-th (every one of the
        # first 100, then every 13th; a call makes its own allocations
        # early, the rest are list items alike)
        for start in [*range(1, 101), *range(101, 20000, 13)]:
            try:
                fail(start, start + 1)
                try:
                    r, x = f(*a), None
                finally:
                    unhook()
            except BaseException as e:
                r, x = None, e
            unhook()
            n += 1
            if isinstance(x, MemoryError):
                continue
            if x is not None and not isinstance(x, (ValueError, TypeError)):
                sys.exit(f"case {i} ({c['fn']}): {name}: allocation {start}: {x!r}")
            got = repr(x) if x is not None else repr(r)
            if got != want:
                sys.exit(f"case {i} ({c['fn']}): {name}: allocation {start}: {got} instead of {want}")
            break
# the table and picture entry points, on a small table
t = {"g": ["b", "a", None], "x": [2.0, float("nan"), 1.0], "s": ["x,y", "q", "z"]}
calls = [("show", (t, 10, 5)), ("show", (t, 1, 1)), ("vshow", ({"a": 1.5, "b": 2.0}, 80)),
         ("csv", (t,)), ("json", (t,)), ("order", ([t["g"], t["x"]], ["g", "x"], [False, True])),
         ("test", (t["x"], "x", ">", 1.0)), ("test", (t["g"], "g", "in", ["a", "b"])),
         ("pairs", ([t["g"]], [["a", "b", "b"]], "full")),
         ("svg", ([([0.0, 1.0], [1.0, 0.5], True, 0, "a"), ([0.5], [0.2], False, 1, "b")],
                  "t", "x", "y", True, 0.0, 2))]
for name, a in calls:
    f = real[name]
    want = repr(f(*a))
    sys.stderr.write(f"table call {name}\\n")
    sys.stderr.flush()
    for start in range(1, 400):
        try:
            fail(start, start + 1)
            try:
                r, x = f(*a), None
            finally:
                unhook()
        except BaseException as e:
            r, x = None, e
        unhook()
        n += 1
        if isinstance(x, MemoryError):
            continue
        if x is not None or repr(r) != want:
            sys.exit(f"table call {name}: allocation {start}: {x!r} {r!r}")
        break
print(n)
"""


@pytest.mark.skipif(os.environ.get("CLAB_NOMEM") != "1",
                    reason="slow (minutes): set CLAB_NOMEM=1 to run it")
def test_every_failed_allocation_is_a_memory_error():
    try:
        import _testcapi
    except ImportError:
        pytest.skip("needs _testcapi")
    if not hasattr(_testcapi, "set_nomemory"):
        pytest.skip("needs _testcapi.set_nomemory")
    r = subprocess.run([sys.executable, "-c", NOMEM], cwd=HERE, capture_output=True,
                       text=True, check=False, timeout=1800)
    assert r.returncode == 0, r.stderr[-3000:]
    assert int(r.stdout) > len(CASES)


class Seq:
    """A sequence whose items cannot be listed."""

    def __len__(self):
        return 2

    def __getitem__(self, i):
        return 1.0

    def __iter__(self):
        raise RuntimeError("no")


class Ai(list):
    """A list that offers an array interface of its own."""

    def __init__(self, v, ai):
        super().__init__(v)
        self.__array_interface__ = ai


def test_odd_sequences_and_array_interfaces():
    with pytest.raises(TypeError, match="`done` must be finite"):
        crimlab.losbnd(Seq(), 1)
    with pytest.raises(TypeError, match="`rel` must be TRUE or FALSE"):
        crimlab.sellab([1, 0], Seq(), [True, True])
    want = crimlab.losbnd([1.0, 2.0], 1)
    # an interface that is no dict, has no address, odd strides or a type
    # string too short is ignored: the list is read
    for ai in (None, {"typestr": "<f8", "shape": (2,), "data": None},
               {"typestr": "<f8", "shape": (2,), "data": (0, True), "strides": (8, 8)},
               {"typestr": "<", "shape": (2,), "data": (0, True)}):
        assert crimlab.losbnd(Ai([1.0, 2.0], ai), 1) == want
    with pytest.raises(MemoryError):
        crimlab.losbnd(Ai([], {"typestr": "<f4", "shape": (2**62,), "data": (8, True)}), 1)
    with pytest.raises(TypeError):
        crimlab.losbnd(Ai([], {"typestr": "<f8", "shape": ("x",), "data": (8, True)}), 1)
    with pytest.raises(TypeError):
        crimlab.losbnd(Ai([], {"typestr": "<f8", "shape": (2,), "data": ("x", True)}), 1)
    with pytest.raises(TypeError):
        crimlab.losbnd(Ai([], {"typestr": "<f8", "shape": (2,), "data": (8, True),
                                "strides": ("x",)}), 1)


def test_a_long_message_is_kept_whole():
    long = "x" * 300
    with pytest.raises(ValueError) as e:
        crimlab.vodtst({"a": 1, "b": 1}, [1, 1], ref=long)
    assert str(e.value) == f"vodtst: `ref` names no group: {long}."


def test_entry_points_check_their_number_of_arguments():
    for name in ("vodtst", "cllr", "selftest", "tlen"):
        with pytest.raises(TypeError, match="expected"):
            getattr(_clab, name)(1, 2, 3, 4, 5, 6, 7)
