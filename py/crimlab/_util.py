# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""What the wrappers do before the compiled call: names and labels.

Numbers are read and checked in the compiled module (ext/clab.cpp). Here a
per-group argument is split into its values and its group names, two named
arguments are aligned by name, and labels (sentences, alleles) are written
as R's as.character() writes them.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping, Sequence
from typing import Any, Protocol

from . import _clab


class Arr(Protocol):
    """Anything NumPy can read as an array (an ndarray, a pandas Series)."""

    def __array__(self, *args: Any, **kwargs: Any) -> Any: ...


Num = int | float
#: Numbers: one number, a sequence (None and NaN are missing), a NumPy
#: array, a pandas Series, or, for a per-group argument, a mapping from
#: group name to number.
Vec = Num | Sequence[Num | None] | Mapping[Any, Num | None] | Arr
#: A data frame: column name to column.
Frame = dict[str, list[Any]]
#: A named numeric vector: name to value.
Named = dict[str, float]

# a treat that matches no label
NEVER = object()


def _pd() -> Any:
    # pandas only when the caller has imported it: it is never a dependency
    return sys.modules.get("pandas")


def series(x: Any) -> bool:
    pd = _pd()
    return pd is not None and isinstance(x, pd.Series)


def dframe(x: Any) -> bool:
    pd = _pd()
    return pd is not None and isinstance(x, pd.DataFrame)


def labelled(idx: Any) -> list[Any] | None:
    """The labels of a pandas index, None for a RangeIndex (positions)."""
    return None if isinstance(idx, _pd().RangeIndex) else list(idx)


def vals(x: Any) -> Any:
    """A Series' values for a numeric argument; others unchanged.

    A categorical is read by its labels, as text: never as numbers.
    """
    if not series(x):
        return x
    if isinstance(x.dtype, _pd().CategoricalDtype):
        return [str(v) for v in x]
    if getattr(x.dtype, "kind", "O") in "fiu":
        return x.to_numpy(dtype="float64", na_value=math.nan)
    return x.to_numpy(dtype=object, na_value=None)


def grp(x: Any) -> tuple[Any, list[Any] | None]:
    """The values and the group names of a per-group argument: the keys of
    a mapping, the index of a Series (unless a RangeIndex), else None."""
    if isinstance(x, Mapping):
        return list(x.values()), list(x.keys())
    if series(x):
        return vals(x), labelled(x.index)
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


def named(fn: str, args: list[tuple[str, Any]],
          strict: bool = True) -> tuple[list[Any], list[Any] | None]:
    """align() on arguments given as they come: a mapping or a Series is
    split into its values and names (grp())."""
    return align(fn, [(nm, *grp(x)) for nm, x in args], strict)


def isbool(v: Any) -> bool:
    return isinstance(v, bool) or type(v).__name__ in ("bool_", "bool")


def lab(v: Any) -> str | None:
    """as.character() of one label: None when missing. TypeError for a
    value that is not text, a number or a bool."""
    if v is None or type(v).__name__ == "NAType":
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
    """A list of x's items: x itself unless it is an array or a Series."""
    if series(x):
        return x.to_numpy(dtype=object, na_value=None).tolist()
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
    (bools among them), "s" otherwise (text, a categorical)."""
    bad = TypeError(f"{fn}: `{nm}` must be labels: a character, numeric or logical vector or a factor.")
    if matrix(z):
        raise notvec(fn, nm)
    cat = series(z) and isinstance(z.dtype, _pd().CategoricalDtype) or \
        type(z).__name__ == "Categorical"
    if hasattr(z, "ravel") and hasattr(z, "__array__") and not series(z):
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
    vals = [v for v in z if v is not None and type(v).__name__ != "NAType"]
    kind = "s"
    if not cat and vals and all(isbool(v) for v in vals):
        kind = "l"
    elif not cat and vals and all(not isinstance(v, (str, bytes)) for v in vals):
        kind = "n"
    return out, kind


def label1(t: Any) -> Any:
    """The label of treat: one text, number or bool; NEVER for anything
    else (or a missing value), which matches no label."""
    if isseq(t) or isinstance(t, Mapping):
        return NEVER
    try:
        v = lab(t)
    except TypeError:
        return NEVER
    return NEVER if v is None else v


def rows(x: Any) -> tuple[list[Any] | None, list[Any] | None]:
    """The rows of a table and its row names: a mapping (row name to row),
    a DataFrame, a 2-D array or a sequence of rows; a sequence of two
    values is one row. Rows are None when x is no table."""
    names = None
    if isinstance(x, Mapping):
        names, r = list(x.keys()), [flat(v) for v in x.values()]
    elif dframe(x):
        names = labelled(x.index)
        r = x.to_numpy(dtype=object, na_value=None).tolist()
    else:
        x = flat(x)
        if not isseq(x):
            return None, None
        r = list(x)
        if len(r) == 2 and not any(isseq(flat(e)) for e in r):
            r = [r]
        else:
            r = [flat(e) for e in r]
    if not all(isseq(e) for e in r):
        return None, names
    return [list(e) for e in r], names


def table(fn: str, nm: str, g: Any) -> tuple[list[str], list[Any] | None]:
    """A genotype table: the alleles row by row as labels, and the loci."""
    r, loci = rows(g)
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
