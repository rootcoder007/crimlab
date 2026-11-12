# crimlab for Python

Criminology and sociolegal statistics: the Python package binds the same
compiled kernels as the R package crimlab (`../src`), with the same
function names, the same argument names and the same results. Every
number in the parity tests (`tests/data/ref.json`, written by
`dev/pyref.R` from the R package) is bit-for-bit what R returns.

There are no runtime dependencies. NumPy arrays and pandas objects are
accepted when you pass them; results are plain Python objects, and
`frame()` / `arrays()` turn them into pandas or NumPy objects when those
are installed. Plots are drawn with matplotlib when it is installed.

## Install

From a clone of the repository (a C++17 compiler is needed):

    pip install ./py

The CI workflow (`.github/workflows/py.yml`) also builds wheels (one per
platform, stable ABI, Python 3.10 and later) and an sdist as workflow
artifacts. The package is not published on PyPI.

## Examples

Disparity:

```python
>>> from crimlab import vodtst, dspbm
>>> r = vodtst(day={"white": 600, "black": 400}, night={"white": 700, "black": 300})
>>> r["group"], [round(x, 4) for x in r["ratio"]]
(['white', 'black'], [1.0, 1.5556])
>>> dspbm(pop={"a": 9000, "b": 1000}, contact={"a": 900, "b": 300},
...       force={"a": 9, "b": 6})["rd"]
[1.0, 6.0]
```

Forensics:

```python
>>> from crimlab import gtmp, cllrmin, kinidx
>>> round(gtmp(p={"D3": 0.1, "vWA": 0.2}, q=[0.3, None], hom=[False, True],
...            theta=0.03)["log10lr"], 4)
2.3068
>>> kinidx({"D3": [12, 14]}, {"D3": [12, 14]}, {"D3": [0.1, 0.2]}, rel="sibling")["locus"]
{'locus': ['D3'], 'idx': [8.375]}
>>> {k: round(v, 4) for k, v in cllrmin([50, 200, 8], [0.01, 0.2, 1.5]).items()}
{'cllr': 0.3008, 'min': 0.0, 'cal': 0.3008}
```

Policing:

```python
>>> from crimlab import clrsmp, expbnd
>>> clrsmp(clra={"minor": 90, "serious": 20}, reca=[100, 50],
...        clrb=[19, 60], recb=[20, 140])["total"]["simpson"]
2.0
>>> expbnd(y={"white": 100, "black": 300}, m=[1000, 1000], gamma=1.5)["sign"]
[0, 1]
```

Corrections:

```python
>>> from crimlab import incap, sntbnd
>>> incap(lam={"low": 2, "high": 10}, q=0.1, s=1, share=[0.8, 0.2])["group"]
['low', 'high', 'all']
>>> sntbnd([1, 0, 1, 1, 0, 0], ["a", "a", "b", "b", "b", "a"], treat="b")["term"]
['a', 'b', 'effect']
```

Proofs and the self-test:

```python
>>> from crimlab import proof, selftest
>>> proof("orrisk")["theorem"][0]
'Crimlab.RelRisk.or_eq_rr_mul'
>>> selftest()["nfail"]
0
```

## Results

| R | Python |
|---|---|
| data frame | dict of columns (lists), in R's column order |
| named numeric vector | dict of floats |
| list | dict |
| `NA` | `None` |
| `NaN`, `Inf` | `float("nan")`, `float("inf")` |

`crimlab.frame(x)` gives a pandas DataFrame (or a Series for a dict of
numbers); `crimlab.arrays(x)` gives NumPy arrays, with `None` as NaN in a
numeric column.

## Inputs

**Numbers** are one number (Python has no vector of one), a sequence
(list, tuple, `array.array`, `memoryview`, `range`), a NumPy array or a
pandas Series. A NumPy array of float64 that is contiguous is read in
place, without a copy (through NumPy's array interface); other numeric
dtypes (float32, any integer type), non-contiguous arrays, byte orders not
the machine's and object arrays are converted. A 2-D array with one row or
one column is the vector it holds; an array with more than one dimension
above 1 is refused with R's message ("must be a vector, not a matrix ...").
A nested list is not flattened: it is refused as not numbers.

Not numbers (TypeError): text, bytes, bools (as in R, a logical is not a
number), dates, complex numbers, dicts and sets where numbers are expected,
iterators and generators, nested sequences.

**Missing values** are `None`, `float("nan")` and `pandas.NA`, all read as
R's `NA`. Where R refuses `NA` they are refused with R's message
(ValueError); where R accepts it (`q` at a homozygous locus of `gtmp()`,
`y` for a detained defendant in `sellab()`), so does Python.

**Logical** arguments (`hom`, `rel`, `rule`) take `True` and `False`
(Python's or NumPy's), one or a sequence or a bool array; numbers are
refused as in R.

**Labels** (the sentences `z` of `sntbnd()`, alleles) are text, numbers or
bools, written as R's `as.character()` writes them: `2.0` is `"2"`,
`0.1 + 0.2` is `"0.3"`, `1e5` is `"1e+05"`, `True` is `"TRUE"`. A
categorical (pandas) is read by its labels, never its integer codes; a
categorical given where numbers are expected is refused, as R refuses a
factor.

**Genotypes** are tables with two columns, one row per locus: a sequence of
pairs, a 2-D array, a DataFrame (its index names the loci unless it is a
RangeIndex) or a mapping from locus to pair; two values are one locus.

## Names

Names are the identity of a group, under the rule of the R package
(`src/bind.cpp`, `align()`). An argument is named when it is a mapping
(dict) or a pandas Series whose index is not a RangeIndex; a DataFrame's
index names its rows.

- All the arguments of a call that hold one value per group (or per unit,
  per locus, per arm) named, with the same names: each is matched to the
  first by name, in any order. `vodtst({"a": 6, "b": 4}, {"b": 5, "a": 5})`
  pairs `a` with `a`.
- All named, with different names: an error that lists them,
  ``vodtst: `night` must name the same groups as `day`; missing: b; extra: c.``
- Some named, some not: positions, as given. Only names protect against a
  mismatch.
- A name given twice in a per-group argument is an error
  (``vodtst: `day` names a group twice: a.``); in a per-unit argument
  (`outtst()`, `lrpost()`, `rtmhs()`, `sellab()`, `sntbnd()`) it leaves
  positions in use.
- Only arguments of the length of the first take part: one value to be
  recycled keeps its place.
- The groups are named by the first argument that names them, else `"1"`,
  `"2"`, ...; the output keeps the order of the first argument. `ref` is a
  position from 1 in that order when it is a number, a group name
  otherwise, and is never chosen by sorting; `vodtst()`, `dspbm()` and
  `expbnd()` mark it in their `ref` column.
- `vodtst(day=table)`: with `night` left out, `day` is one table whose two
  columns are named `day` and `night`, a DataFrame (its index names the
  groups) or a dict of the two columns.
- `dyadrt()`: a DataFrame `count` that labels its rows and its columns has
  its columns matched to its rows by name, and both to `pop` when `pop` is
  named.
- Genotype tables and `freq` that all name their loci are matched by locus
  to the first.
- `hzdsel()`: `surh` and `surl` both named `control` and `treated`, in
  either order, are read by name; otherwise the first value is the control
  arm.

## Errors

Errors are exceptions whose message is R's, word for word, starting with
the function's name: `ValueError` for a value outside what the function
takes (missing, not finite, out of range, wrong length), `TypeError` for an
argument of the wrong kind. A failed allocation is a `MemoryError`.

## Differences from R

1. `orrisk()`'s first argument is `or`, as in R, but `or` is a Python
   keyword: give it by position, or as `orrisk(**{"or": 3}, rate=0.4,
   share=0.5)`.
2. Results are dicts and lists (see Results). An unnamed result of one
   number is a float for `cllr()`; `lrpost()` always returns a list.
3. R's forms of a named vector (a 1-d table, a list of single numbers, a
   one-row data frame) are a dict or a Series in Python.
4. `ref`: a bool or None is refused with a TypeError; a number is always a
   position, so a group named by a number cannot be chosen by name. R
   reads a factor `ref` by its label; Python has no one-value categorical.
5. `tipplt()` and `eceplt()` return the plotted data (R returns it
   invisibly) and draw only when matplotlib is installed; `draw=False`
   skips drawing, `ax` picks the axes and other keyword arguments go to
   `ax.set()` (R passes `...` to `plot()`).
6. None, NaN and pandas.NA are all R's `NA` on input. On output R's `NA`
   is None and `NaN` is NaN.
7. A nested list where a vector is expected is refused as not numbers
   (TypeError); R's message for a matrix is given for NumPy arrays with
   more than one dimension above 1. For the allele names of `mixlr()` a
   2-D array is refused as not names of `p`.
8. `treat` in `sntbnd()` must be one label; a list, even of one label, is
   refused (R accepts a vector of length 1). The default (True for bools,
   1 for numbers 0 and 1) follows R.
9. A one-row genotype table is named `"1"` unless it is named (mapping
   key, DataFrame index); R names a one-row matrix by its only column name
   when it has column names but no row names.
10. `proof()`: `None` means every function (R's `NULL`); R's `NA` has no
    Python spelling, and `[None]` is refused. `proof()` lists the network
    functions too, which are not bound in Python yet.
11. When several arguments are wrong at once the first error reported may
    differ: names and labels are checked in Python before the numbers are
    checked in the compiled module.
12. R conditions carry a call; Python exceptions have none.
13. `selftest()` is exported (R runs the self-test through an internal
    `.Call`).
14. R-only kinds of input have no counterpart: `"bytes"`-encoded strings,
    classed numbers with units, labelled vectors (R refuses them with a
    message of its own), environments.
15. Not bound yet: the network functions `ntdeg`, `ntbtw`, `ntcls`,
    `ntevc`, `ntcore`, `ntcmp`, `ntlvn`, `ntmod`, `ntbrk`, `ntchg` and
    `coofg`. Their kernels are compiled into the module and run in
    `selftest()`.

## Development

From the repository root:

    pip install ./py pytest
    cd py/tests && python -m pytest
    python -m pytest --doctest-modules --pyargs crimlab
    Rscript dev/pyref.R LIB     # regenerate tests/data/ref.json from R

`setup.py` copies the kernel sources (`../src/*.cpp` and `*.h` that do not
include R), the version (`../DESCRIPTION`), the proof table
(`../src/bind.cpp`), the Lean files (`../inst/lean`) and the licence into
`py/` before every build; those copies are build output and are not
committed. An sdist carries them and builds on its own.

Adding a function: its entry point in `ext/clab.cpp` (one function and one
line in the method table), its wrapper with a docstring in `crimlab/`, its
name in `__all__` and `_clab.pyi`, and its cases in `dev/pyref.R`.

A documentation site is left to the maintainer; the docstrings follow the
numpydoc layout, so Sphinx with `numpydoc` or `napoleon` can render them.
