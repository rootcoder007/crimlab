# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Results: dicts that print like R and export themselves.

A table is a dict of equal-length columns (lists), a vector a dict of
numbers, a record a dict of results; each is still the dict it was, so
``r["ratio"]`` and ``dict(r)`` work as before. Printing, CSV, JSON, sort,
filter and join run in the compiled module.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from typing import Any

from . import _clab


class Table(dict[str, list[Any]]):
    """A data frame: a dict of columns, all lists of one length.

    Printing gives an aligned table (row numbers first, numbers to 7
    significant digits, a long table cut to its head and tail with the row
    count). ``Table(d)`` makes one from any dict of equal-length columns.
    Beyond the dict: ``shape``, ``rows()``, ``row(i)``, ``csv()``,
    ``json()``, ``sort()``, ``filter()`` and ``join()``, each returning a
    new table or text and leaving this one unchanged.

    Examples
    --------
    >>> t = Table({"g": ["a", "b", "c"], "x": [3.0, 1.5, None]})
    >>> t
      g   x
    1 a 3.0
    2 b 1.5
    3 c  NA
    >>> t.sort("x").row(0), t.shape
    ({'g': 'b', 'x': 1.5}, (3, 2))
    """

    def __init__(self, cols: Mapping[str, Iterable[Any]] | None = None, /) -> None:
        super().__init__((k, list(v)) for k, v in (cols or {}).items())
        if len({len(v) for v in self.values()}) > 1:
            raise ValueError("Table: the columns must have one length.")

    def __repr__(self) -> str:
        return str(_clab.show(self, 10, 5))

    __str__ = __repr__

    @property
    def shape(self) -> tuple[int, int]:
        """(rows, columns)."""
        return (len(next(iter(self.values()))) if self else 0, len(self))

    def rows(self) -> list[dict[str, Any]]:
        """The rows, each a dict of column name to value."""
        return [dict(zip(self, r)) for r in zip(*self.values())]

    def row(self, i: int) -> dict[str, Any]:
        """Row i (from 0; negative from the end) as a dict."""
        return {k: v[i] for k, v in self.items()}

    def csv(self) -> str:
        """CSV text (RFC 4180): a header row, CRLF line ends, fields with a
        comma, a quote or a line break quoted. Numbers in full (the shortest
        text that reads back to the same double), None (R's NA) as an
        empty field, NaN and the infinities as NaN, Inf and -Inf, bools as
        TRUE and FALSE."""
        return str(_clab.csv(self))

    def json(self) -> str:
        """JSON text: an object of column arrays. None (R's NA) is null;
        NaN and the infinities, which JSON has no numbers for, are the
        strings "NaN", "Inf" and "-Inf"."""
        return str(_clab.json(self))

    def sort(self, by: str | list[str], desc: bool | list[bool] = False) -> Table:
        """The rows in order of one column or several (a list, the first
        deciding first). Stable; ``desc`` one bool or one per column; None
        and NaN last in either direction; text by Unicode code point,
        numbers (bools as 0 and 1) numerically; a column of text and
        numbers both is an error."""
        cols = [by] if isinstance(by, str) else list(by)
        ds = [bool(desc)] * len(cols) if isinstance(desc, bool) else list(desc)
        if len(ds) != len(cols):
            raise ValueError("sort: `desc` must be one bool or one per column of `by`.")
        for c in cols:
            if c not in self:
                raise ValueError(f"sort: no column `{c}`.")
        o = _clab.order([self[c] for c in cols], cols, ds)
        return Table({k: [v[i] for i in o] for k, v in self.items()})

    def filter(self, f: Callable[[dict[str, Any]], Any] | None = None, /,
               **conds: Any) -> Table:
        """The rows where every condition holds. A keyword condition names
        a column: ``col=(op, value)`` with op one of ``==``, ``!=``, ``<``,
        ``<=``, ``>``, ``>=`` and ``in`` (value a list), or ``col=value``
        for ``==``. f, a function of a row (a dict), keeps the rows where it
        is true. A missing value (None, NaN) never passes; text and numbers
        are not ordered against each other (an error)."""
        keep = [True] * self.shape[0]
        for c, v in conds.items():
            if c not in self:
                raise ValueError(f"filter: no column `{c}`.")
            op, val = v if isinstance(v, tuple) and len(v) == 2 and isinstance(v[0], str) \
                else ("==", v)
            keep = [a and b for a, b in zip(keep, _clab.test(self[c], c, op, val))]
        if f is not None:
            if not callable(f):
                raise TypeError("filter: the first argument must be a function of a row.")
            keep = [a and bool(f(r)) for a, r in zip(keep, self.rows())]
        return Table({k: [x for x, y in zip(v, keep) if y] for k, v in self.items()})

    def join(self, other: Mapping[str, Iterable[Any]], on: str | list[str],
             how: str = "inner") -> Table:
        """The rows of this table and another (any dict of equal-length
        columns) matched on the key columns ``on``. ``how``: inner (pairs
        only), left (and the rows of this table without a match), right
        (and the other's), full (both). Keys match by exact equality, text
        with text, numbers with numbers, bools with bools; a missing key
        (None, NaN) matches nothing. Rows: this table's order, each with its
        matches in the other's order, then the other's unmatched rows (right,
        full); a key on n rows here and m there gives n * m rows. Columns:
        the keys once, then this table's, then the other's; a name in both
        that is not a key gets ".x" (here) and ".y" (there). Unmatched
        cells are None."""
        o = other if isinstance(other, Table) else Table(other)
        keys = [on] if isinstance(on, str) else list(on)
        for k in keys:
            for t, side in ((self, "this"), (o, "the other")):
                if k not in t:
                    raise ValueError(f"join: no key column `{k}` in {side} table.")
        li, ri = _clab.pairs([self[k] for k in keys], [o[k] for k in keys], how)
        out: dict[str, list[Any]] = {}
        for k in keys:
            out[k] = [self[k][i] if i >= 0 else o[k][j] for i, j in zip(li, ri)]
        for k, v in self.items():
            if k not in keys:
                out[k + ".x" if k in o else k] = [v[i] if i >= 0 else None for i in li]
        for k, v in o.items():
            if k not in keys:
                out[k + ".y" if k in self else k] = [v[j] if j >= 0 else None for j in ri]
        return Table(out)


class Vector(dict[str, Any]):
    """A named vector: a dict of name to number. Prints as R prints a
    named vector, names above values.

    Examples
    --------
    >>> Vector({"stop": 0.5, "hit": 0.75})
    stop  hit
    0.50 0.75
    """

    def __repr__(self) -> str:
        return str(_clab.vshow(self, 80))

    __str__ = __repr__

    def json(self) -> str:
        """JSON text: an object, by the rule of :meth:`Table.json`."""
        return str(_clab.json(self))


class Record(dict[str, Any]):
    """A list of results: a dict of name to table, vector or number.
    Prints each part under ``$name``, as R prints a list."""

    def __repr__(self) -> str:
        def one(v: Any) -> str:
            if isinstance(v, float) and not isinstance(v, bool):
                return "[1] " + str(_clab.vshow({"": v}, 80)).split("\n")[-1].strip()
            return repr(v)

        return "\n\n".join(f"${k}\n{one(v)}" for k, v in self.items())

    __str__ = __repr__

    def json(self) -> str:
        """JSON text: an object, by the rule of :meth:`Table.json`."""
        return str(_clab.json(self))


def tab(x: Mapping[str, Any]) -> Table:
    return Table(x)


def vec(x: Mapping[str, Any]) -> Vector:
    return Vector(x)


def rec(x: Mapping[str, Any]) -> Record:
    return Record({k: wrap(v) for k, v in x.items()})


def wrap(x: Any) -> Any:
    """A result of the compiled module as a Table, a Vector or a Record."""
    if not isinstance(x, dict):
        return x
    v = list(x.values())
    if v and all(isinstance(c, list) for c in v) and len({len(c) for c in v}) == 1:
        return Table(x)
    if all(c is None or isinstance(c, (int, float)) for c in v):
        return Vector(x)
    return Record({k: wrap(c) for k, c in x.items()})
