# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""What the wrappers do before the compiled call: names and labels.

Numbers are read and checked in the compiled module (ext/clab.cpp). Here a
per-group argument is split into its values and its group names, two named
arguments are aligned by name, and labels (sentences, alleles) are written
as R's as.character() writes them.
"""

from __future__ import annotations

import json
import keyword
import math
from collections.abc import Mapping, Sequence
from typing import Any, Protocol

from . import _clab


class Arr(Protocol):
    """An array: an object that offers __array_interface__ (or __array__)."""

    def __array__(self, *args: Any, **kwargs: Any) -> Any: ...


Num = int | float
#: Numbers: one number, a sequence (None and NaN are missing), an array,
#: or, for a per-group argument, a mapping from group name to number or a
#: sequence of (name, number) tuples (a name may repeat there).
Vec = Num | Sequence[Num | None] | Mapping[Any, Num | None] | Sequence[tuple[Any, Num | None]] | Arr
#: A data frame: column name to column.
Frame = dict[str, list[Any]]
#: A named numeric vector: name to value.
Named = dict[str, float]

# a treat that matches no label
NEVER = object()


def items(x: Any) -> list[tuple[Any, Any]] | None:
    """The (name, value) pairs of a mapping, in iteration order: any object
    but text that has a callable items() (a dict, any Mapping, an object
    such as a pandas Series that yields (label, value)); None otherwise.
    Names may repeat."""
    if isinstance(x, (str, bytes)):
        return None
    f = getattr(x, "items", None)
    return list(f()) if callable(f) else None


def mapping(x: Any) -> bool:
    return not isinstance(x, (str, bytes)) and callable(getattr(x, "items", None))


def scalar(v: Any) -> bool:
    return not isseq(v) and not mapping(v) and not hasattr(v, "__array__")


def pairs(x: Any) -> bool:
    """True for a list or tuple of (name, value) tuples: named values whose
    names may repeat, which a dict cannot hold."""
    return isinstance(x, (list, tuple)) and len(x) > 0 and all(
        isinstance(e, tuple) and len(e) == 2 and scalar(e[0]) and scalar(e[1]) for e in x)


def grp(x: Any) -> tuple[Any, list[Any] | None]:
    """The values and the group names of a per-group argument: the pairs
    of a mapping (items()), or of a sequence of (name, value) tuples; no
    names otherwise."""
    p = items(x)
    if p is None and pairs(x):
        p = list(x)
    if p is not None:
        return [v for _, v in p], [k for k, _ in p]
    return x, None


def size(x: Any) -> int:
    """The number of values of an argument; 1 for a single one."""
    try:
        return len(x)
    except TypeError:
        return 1


def twice(n: list[Any]) -> Any:
    """The first name that repeats in n, else None."""
    seen: set[Any] = set()
    for k in n:
        if k in seen:
            return k
        seen.add(k)
    return None


def order(fn: str, an: str, a: list[Any], bn: str, b: list[Any], v: list[Any]) -> list[Any]:
    """The values v of argument bn (names b) in the order of the names a of
    argument an; the two sets of names must be the same."""
    sa, sb = set(a), set(b)
    if sa != sb:
        miss = ", ".join(str(k) for k in a if k not in sb)
        extra = ", ".join(str(k) for k in b if k not in sa)
        raise ValueError(f"{fn}: {bn} must name the same groups as {an}; missing: {miss}; "
                         f"extra: {extra}.")
    at = {k: i for i, k in enumerate(b)}
    return [v[at[k]] for k in a]


def align(fn: str, args: list[tuple[str, Any, list[Any] | None]],
          strict: bool) -> tuple[list[Any], list[Any] | None]:
    """R's rule for the arguments of one call that hold one value per group
    (src/bind.cpp, align()). args holds (name, values, names) for each.
    Only arguments of the length of the first take part. All of them named
    with one set of names: each later one is put in the order of the first.
    Named with different sets: an error listing the names missing and extra.
    Some named, some not: positions. A name given twice is an error (strict)
    or leaves positions in use (not strict, for per-unit arguments). The
    groups are named by the first named argument, else None."""
    vs = [v for _, v, _ in args]
    n = size(vs[0])
    use = [v is not None and size(v) == n for v in vs]
    named = [j for j, (_, _, nm) in enumerate(args) if use[j] and nm is not None]
    dup = False
    for j in named:
        k = twice(args[j][2] or [])
        if k is not None:
            if strict:
                raise ValueError(f"{fn}: `{args[j][0]}` names a group twice: {k}.")
            dup = True
    if len(named) == sum(use) > 1 and not dup:
        a0 = args[0]
        for j in range(1, len(args)):
            if use[j]:
                vs[j] = order(fn, f"`{a0[0]}`", a0[2] or [], f"`{args[j][0]}`",
                              args[j][2] or [], list(vs[j]))
    return vs, args[named[0]][2] if named else None


def byname(fn: str, args: list[tuple[str, Any]]) -> list[tuple[Any, list[Any] | None]]:
    """R's rule for arguments with one value per unit (src/spin.cpp,
    byname()): when every argument given (not None) is named, with names
    that do not repeat, all of one length, and at least two are given, each
    later one is put in the order of the first, which must name the same
    units; otherwise positions are used. Returns (values, names) for each."""
    vs = [grp(x) if x is not None else (None, None) for _, x in args]
    given = [j for j, (v, _) in enumerate(vs) if v is not None]
    if len(given) < 2 or any(vs[j][1] is None for j in given):
        return vs
    n = size(vs[given[0]][0])
    if any(size(vs[j][0]) != n or twice(vs[j][1] or []) is not None for j in given):
        return vs
    f = given[0]
    fn_ = vs[f][1] or []
    for j in given[1:]:
        at = {k: i for i, k in enumerate(vs[j][1] or [])}
        for k in fn_:
            if k not in at:
                raise ValueError(f"{fn}: `{args[j][0]}` must name the same units as "
                                 f"`{args[f][0]}`; `{args[j][0]}` has no {k}.")
        vv = list(vs[j][0])
        vs[j] = ([vv[at[k]] for k in fn_], list(fn_))
    return vs


def named(fn: str, args: list[tuple[str, Any]],
          strict: bool = True) -> tuple[list[Any], list[Any] | None]:
    """align() on arguments given as they come: a mapping or a sequence of
    (name, value) tuples is split into its values and names (grp())."""
    return align(fn, [(nm, *grp(x)) for nm, x in args], strict)


def isbool(v: Any) -> bool:
    return isinstance(v, bool) or type(v).__name__ in ("bool_", "bool")


def lab(v: Any) -> str | None:
    """as.character() of one label: None when missing. TypeError for a
    value that is not text, a number or a bool."""
    if v is None:
        return None
    if isinstance(v, str):
        return str(v)
    if isbool(v):
        return "TRUE" if v else "FALSE"
    if hasattr(v, "__index__"):
        return str(v.__index__())
    if isinstance(v, float) or type(v).__name__.startswith("float"):
        f = float(v)
        return None if math.isnan(f) else str(_clab.rchar(f))
    raise TypeError(type(v).__name__)


def flat(x: Any) -> Any:
    """A list of x's items: x itself unless it is an array (.tolist())."""
    return x.tolist() if hasattr(x, "tolist") and hasattr(x, "__array__") else x


def isseq(x: Any) -> bool:
    return isinstance(x, Sequence) and not isinstance(x, (str, bytes))


def matrix(x: Any) -> bool:
    """True for an array with more than one dimension above 1."""
    return sum(1 for d in getattr(x, "shape", ()) if d > 1) > 1


def notvec(fn: str, nm: str) -> ValueError:
    return ValueError(f"{fn}: `{nm}` must be a vector, not a matrix with more than one row and column.")


def labels(fn: str, nm: str, z: Any) -> tuple[list[str | None], str]:
    """The labels of z, as.character() of each (None where missing), and
    its kind: "l" when every label is a bool, "n" when they are numbers
    (bools among them), "s" otherwise."""
    bad = TypeError(f"{fn}: `{nm}` must be labels: a character, numeric or logical vector or a factor.")
    if matrix(z):
        raise notvec(fn, nm)
    if hasattr(z, "ravel") and hasattr(z, "__array__"):
        z = z.ravel()
    z = flat(z)
    if isinstance(z, (str, bytes)):
        z = [z]
    if not isseq(z):
        raise bad
    try:
        out = [lab(v) for v in z]
    except TypeError:
        raise bad from None
    vals = [v for v in z if v is not None]
    kind = "s"
    if vals and all(isbool(v) for v in vals):
        kind = "l"
    elif vals and all(not isinstance(v, (str, bytes)) for v in vals):
        kind = "n"
    return out, kind


def label1(t: Any) -> Any:
    """The label of treat: one text, number or bool; NEVER for anything
    else (or a missing value), which matches no label."""
    if isseq(t) or mapping(t):
        return NEVER
    try:
        v = lab(t)
    except TypeError:
        return NEVER
    return NEVER if v is None else v


def row(v: Any, keys: list[Any] | None) -> Any:
    """One row of a table: a sequence, or a mapping read in the order keys
    (None when its keys are not those)."""
    p = items(v)
    if p is not None:
        d = dict(p)
        return [d[k] for k in keys] if keys is not None and set(d) == set(keys) else None
    return flat(v)


def bycol(p: list[tuple[Any, Any]]) -> tuple[list[Any] | None, list[Any] | None, list[Any]]:
    """A table given as a mapping of columns, each a mapping from row label
    to value (a pandas DataFrame reads so through items()): its rows (None
    when the columns do not label the same rows), row names and columns."""
    cols = [k for k, _ in p]
    lab = [items(c) or [] for _, c in p]
    names = [k for k, _ in lab[0]] if lab else []
    out = []
    for c in lab:
        ks = [k for k, _ in c]
        if ks == names:
            out.append([v for _, v in c])
        elif sorted(map(repr, ks)) == sorted(map(repr, names)) and len(set(ks)) == len(ks):
            d = dict(c)
            out.append([d[k] for k in names])
        else:
            return None, names, cols
    return [list(r) for r in zip(*out)] if out else [], names, cols


def rows(x: Any) -> tuple[list[Any] | None, list[Any] | None, list[Any] | None]:
    """The rows of a table, its row names and its column names. A table is
    a mapping of columns, each a mapping from row label to value (row and
    column names); a mapping from row name to row; a sequence of (row name,
    row) tuples (row names that repeat); a 2-D array; or a sequence of rows,
    each a sequence or a mapping from column name to value (every row
    naming the same columns). A sequence of two values is one row. Rows are
    None when x is no table."""
    names = None
    p = items(x)
    if p is not None and p and all(mapping(v) for _, v in p):
        return bycol(p)
    if p is not None:
        names, r = [k for k, _ in p], [v for _, v in p]
    elif isinstance(x, (list, tuple)) and x and all(
            isinstance(e, tuple) and len(e) == 2 and scalar(e[0]) and not scalar(e[1])
            for e in x):
        names, r = [e[0] for e in x], [e[1] for e in x]
    else:
        x = flat(x)
        if not isseq(x):
            return None, None, None
        r = list(x)
        if len(r) == 2 and not any(isseq(flat(e)) or mapping(e) for e in r):
            r = [r]
    first = items(r[0]) if r else None
    cols = [k for k, _ in first] if first is not None else None
    r = [row(e, cols) for e in r]
    if not all(isseq(e) for e in r):
        return None, names, cols
    return [list(e) for e in r], names, cols


def table(fn: str, nm: str, g: Any) -> tuple[list[str], list[Any] | None]:
    """A genotype table: the alleles row by row as labels, and the loci."""
    r, loci, _ = rows(g)
    msg = f"{fn}: `{nm}` must be a two-column matrix of alleles, one row per locus, with no missing values."
    if r is None or any(len(e) != 2 for e in r):
        raise ValueError(msg)
    out = []
    for v in (v for e in r for v in e):
        if isbool(v):
            raise TypeError(msg)
        try:
            a = lab(v)
        except TypeError:
            raise TypeError(msg) from None
        if a is None:
            raise ValueError(msg)
        out.append(a)
    return out, loci


def nrows(fn: str, nm: str, v: list[Any], n: int) -> None:
    if len(v) // 2 != n:
        raise ValueError(f"{fn}: `{nm}` must have {n} rows, not {len(v) // 2}.")


def jarg(x: Any) -> Any:
    """An argument as the JSON value clab_jrun() reads: a mapping (anything
    with items()) an object, a sequence or an array an array, NaN null (R's
    NA), an infinity the text R would refuse."""
    if x is None or isinstance(x, (bool, str)):
        return x
    if isinstance(x, (int, float)):
        if isinstance(x, float) and math.isnan(x):
            return None
        if isinstance(x, float) and math.isinf(x):
            return "Inf" if x > 0 else "-Inf"
        return x
    p = items(x)
    if p is not None:
        return {str(k): jarg(v) for k, v in p}
    if hasattr(x, "tolist"):
        return jarg(x.tolist())
    try:
        return [jarg(v) for v in x]
    except TypeError:
        raise TypeError(f"cannot pass a {type(x).__name__} to crimlab") from None


def jnum(x: Any) -> Any:
    """NaN and the infinities, written as text in JSON, as floats."""
    if isinstance(x, str) and x in ("NaN", "Inf", "-Inf"):
        return float(x)
    if isinstance(x, list):
        return [jnum(v) for v in x]
    if isinstance(x, dict):
        return {k: jnum(v) for k, v in x.items()}
    return x


def jcall(fn: str, **kw: Any) -> Any:
    """fn of the R package through clab_jrun(), with the arguments given
    (None: left out); the result as Python objects (crimlab._out.wrap)."""
    from ._out import wrap

    # an R name that is a Python keyword (from, lambda) has a final underscore
    a = {(k[:-1] if k.endswith("_") and keyword.iskeyword(k[:-1]) else k): jarg(v)
         for k, v in kw.items() if v is not None}
    r = jnum(json.loads(_clab.jrun(fn, json.dumps(a, allow_nan=False)), parse_int=float))
    return wrap(r)
