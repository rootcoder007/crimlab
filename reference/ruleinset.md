# Ready-made rules

The checks that account for most of what a table is asked: a value in a
set, a number in a range, no missing value, no repeated identifier, a
pattern, an order, an outlier, complete rows, distinct rows, a number of
columns. Each returns a rule for
[`valrules()`](https://rootcoder007.github.io/crimlab/reference/rul.md),
as [`rul()`](https://rootcoder007.github.io/crimlab/reference/rul.md)
does.

## Usage

``` r
ruleinset(col, set, napass = TRUE, sev = c("warning", "fatal"))

rulebtw(col, lo, hi, napass = TRUE, sev = c("warning", "fatal"))

rulenull(col, sev = c("warning", "fatal"))

ruleuniq(col, sev = c("warning", "fatal"))

rulerx(col, pat, napass = TRUE, sev = c("warning", "fatal"))

ruleincr(col, strict = FALSE, sev = c("warning", "fatal"))

rulemad(col, n = 3, napass = TRUE, sev = c("warning", "fatal"))

rulecmpl(sev = c("warning", "fatal"))

ruledstn(cols = NULL, sev = c("warning", "fatal"))

rulecolcnt(n, sev = c("warning", "fatal"))
```

## Arguments

- col:

  The column the rule reads: one non-empty string.

- set:

  The allowed values.

- napass:

  `TRUE` (the default) lets a missing value pass.

- sev:

  `"warning"` (the default) or `"fatal"`: how serious a failure is.
  [`valrules()`](https://rootcoder007.github.io/crimlab/reference/rul.md)
  records it; the caller decides what follows.

- lo, hi:

  The bounds, both included.

- pat:

  A regular expression.

- strict:

  `TRUE` for strictly increasing.

- n:

  `rulemad()`: the number of median absolute deviations allowed.
  `rulecolcnt()`: the number of columns.

- cols:

  The columns over which rows must be distinct, or `NULL` (all of them).

## Value

A rule: an object of class `clrule`, a list of `rule` (the kind),
`name`, `col` (for a rule on one column), `sev`, `msg` and what the rule
was given.

## Details

`ruleinset()`: every value is one of `set` (compared as text).
`rulebtw()`: every value lies in `[lo, hi]`. `rulenull()`: no value is
missing. `ruleuniq()`: no value repeats; the rows reported are the
repeats, not the first occurrence. `rulerx()`: every value matches the
regular expression `pat`. `ruleincr()`: the column never falls (or, with
`strict`, always rises) down the rows; missing values are passed over.
`rulemad()`: every value lies within `n` median absolute deviations of
the median, the deviation scaled by 1.4826 as
[`stats::mad()`](https://rdrr.io/r/stats/mad.html) scales it (so that it
estimates the standard deviation of normal data; Leys et al. 2013 argue
for this screen over one built on the mean and the standard deviation,
which the outliers themselves inflate); when the deviation is 0 (more
than half the values equal) any other value fails. `rulecmpl()`: no row
has a missing value. `ruledstn()`: no row repeats, over all columns or
over `cols`. `rulecolcnt()`: the table has `n` columns.

With `napass = TRUE`, the default, a missing value passes a rule about
values (missingness is `rulenull()`'s question); with `FALSE` it fails.
`rulebtw()`, `ruleincr()` and `rulemad()` read numbers (or dates, or
`TRUE` and `FALSE`): on a column of text they are reported as rules that
could not be evaluated, not passed.

A pattern is read by [`grepl()`](https://rdrr.io/r/base/grep.html) with
`perl = TRUE`; the other languages of the package read it as ECMAScript,
and the syntax the two share (classes, anchors, counts, alternation,
groups) means the same in both. A pattern is at most 256 bytes, and a
value it is tried on at most 4096.

Errors: `col` not one non-empty string; `lo` or `hi` not one number, or
`lo` above `hi`; `n` not one positive number (`rulecolcnt()`: one whole
number of at least 0); `pat` not one string of at most 256 bytes, or not
a regular expression; `cols` neither `NULL` nor column names; `sev` not
`"warning"` or `"fatal"`.

## References

Leys, C., Ley, C., Klein, O., Bernard, P. and Licata, L. (2013).
Detecting outliers: do not use standard deviation around the mean, use
absolute deviation around the median. *Journal of Experimental Social
Psychology*, 49(4), 764-766.
[doi:10.1016/j.jesp.2013.03.013](https://doi.org/10.1016/j.jesp.2013.03.013)

## See also

[`rul()`](https://rootcoder007.github.io/crimlab/reference/rul.md) and
[`valrules()`](https://rootcoder007.github.io/crimlab/reference/rul.md).

## Examples

``` r
d <- data.frame(year = c(2019, 2020, 2021, 2021),
                n = c(410, 395, 9000, 402),
                div = c("11", "14", "14", "5X"))
valrules(d, list(
  rulebtw("year", 2000, 2026),
  ruleincr("year", strict = TRUE),
  rulemad("n", 5),
  rulerx("div", "^[0-9]{2}$"),
  ruledstn(c("year", "div")),
  rulecolcnt(3),
  rulecmpl()
))
#> valrules: 7 rules, 3 failed (0 fatal)
#>   warning  year_increasing: Column 'year' is not strictly increasing
#>   warning  n_within_mads: Column 'n' has values more than 5 MADs from its median (1 of 4 rows); rows 3
#>   warning  div_regex: Column 'div' has values not matching /^[0-9]{2}$/ (1 of 4 rows); rows 4
try(rulebtw("year", 2026, 2000))
#> Error : rulebtw: `lo` must not be above `hi`.
```
