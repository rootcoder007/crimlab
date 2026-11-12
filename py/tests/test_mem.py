# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Memory and references: every call, valid or not, releases what it took,
and a failed allocation anywhere is a MemoryError, never a crash, a wrong
result or a leak. The failing allocations need CPython's _testcapi and
CLAB_NOMEM=1 (they take minutes); the reference counts a debug build
(sys.gettotalrefcount); each test is skipped where its tool is missing."""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from ref import CASES, args

import crimlab
from crimlab import _clab

HERE = Path(__file__).parent


def run_all(k):
    for _ in range(k):
        for c in CASES:
            try:
                getattr(crimlab, c["fn"])(**args(c))
            except (ValueError, TypeError):
                pass


@pytest.mark.skipif(not hasattr(sys, "gettotalrefcount"), reason="needs a debug build")
def test_reference_counts_do_not_grow():
    run_all(2)
    a = sys.gettotalrefcount()
    run_all(10)
    b = sys.gettotalrefcount()
    assert b - a < 100, (a, b)


NOMEM = """
import json, sys, _testcapi, crimlab
from ref import CASES, args
n = 0
# the small cases reach every allocation site; a large one would fail each
# of its thousands of allocations in turn, at a cost quadratic in its size
for c in [c for c in CASES if len(json.dumps(c["args"])) < 800]:
    f = getattr(crimlab, c["fn"])
    a = args(c)
    try:
        want = repr(f(**a))
    except (ValueError, TypeError) as e:
        want = repr(e)
    # every one of the first 100 allocations, then every 13th: a call makes
    # its own allocations early, the rest are list items alike
    for start in [*range(1, 101), *range(101, 20000, 13)]:
        _testcapi.set_nomemory(start, 0)
        try:
            r, x = f(**a), None
        except BaseException as e:
            r, x = None, e
        _testcapi.remove_mem_hooks()
        n += 1
        if isinstance(x, MemoryError):
            continue
        if x is not None and not isinstance(x, (ValueError, TypeError)):
            sys.exit(f"{c['fn']}: allocation {start}: {x!r}")
        got = repr(x) if x is not None else repr(r)
        if got != want:
            sys.exit(f"{c['fn']}: allocation {start}: {got} instead of {want}")
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
    assert r.returncode == 0, r.stderr[-2000:]
    assert int(r.stdout) > 1000


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
