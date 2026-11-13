# Which labelling reproduces a reported contrast?

A reported odds ratio, risk ratio or risk difference depends on which
group is the reference, which outcome is the event and which margin of
the table holds the groups. `relchk()` recomputes the contrast from the
table under every such reading and says which ones reproduce the
reported value: the labelling the report states, another one, or none.

## Usage

``` r
relchk(
  x,
  value,
  type = "or",
  group = NULL,
  ref = NULL,
  event = NULL,
  digits = NULL,
  tol = NULL
)
```

## Arguments

- x:

  The table of counts, groups in rows and outcomes in columns: a matrix,
  a 2-d table, a data frame of numeric columns, or a list of numeric
  vectors of one length (the columns). At least 2 x 2; counts finite and
  at least 0 with no NA, no row or column all 0, and no label missing or
  given twice. Rows and columns without names are labelled "1", "2", ...
  .

- value:

  The reported value: one finite number.

- type:

  The kind of contrast the report claims: `"or"` (odds ratio), `"rr"`
  (risk ratio) or `"rd"` (risk difference).

- group, ref, event:

  The labelling the report claims, each `NULL` (not stated) or a label
  (matched exactly, case included; a factor by its label) or a position:
  `group` and `ref` name rows of `x`, `event` a column.

- digits:

  The decimals `value` was rounded to: a reading matches when it lies
  within half a unit of the last decimal, at either end, so that both
  rounding conventions at the half are accepted. By default the fewest
  decimals that print `value` exactly (give `digits = 2` for a value
  reported as 2.30). Give `digits` or `tol`, not both.

- tol:

  Instead of `digits`, the largest difference allowed: one number
  from 0. Either way, a reading also matches when it differs only in the
  last few bits.

## Value

An object of class `relchk`, a list: `verdict` (`"holds"`: a reading the
claim allows matches; `"other"`: only readings it does not allow match;
`"none"`), `ok` (`TRUE` when the claim holds), `value`, `tol` (the
difference allowed), `digits` (`NA` when `tol` was given or no decimals
print `value` exactly; then `tol` is 0), `type`, `claim` (the stated
`group`, `ref` and `event`, `NA` where not stated) and `reads`, a data
frame with one row per reading: `grp` (`"rows"` or `"cols"`: the margin
holding the groups), `type`, `group`, `ref`, `event`, `value`, `match`
and `claim` (`TRUE` for the readings the claim allows: groups in rows,
the kind `type`, and the stated labels).

## Details

A reading compares group `group` with group `ref` for outcome `event`,
the other outcomes pooled as non-events: with `a` and `b` the events and
non-events of `group`, `c` and `d` those of `ref`, the odds ratio is
`a d / (b c)`, the risk ratio `(a / (a + b)) / (c / (c + d))` and the
risk difference `a / (a + b) - c / (c + d)`. Every ordered pair of
groups and every outcome is read, with the groups in the rows and, the
table transposed, in the columns, and all three kinds of contrast.

No correction is made for zero cells: a reading whose denominator is 0
is `Inf` (a positive numerator) or `NaN` (0 over 0) and matches nothing.
Counts need not be whole numbers (weighted totals are read as they are).

For a 2 x 2 table the odds ratio inverts when the reference group or the
event changes and is unchanged when both change or the table is
transposed; the risk ratio inverts with the reference group, and the
risk difference changes sign with either; `proof("relchk")` names the
Lean theorems.

## Examples

``` r
tab <- matrix(c(20, 10, 80, 90), 2,
              dimnames = list(c("A", "B"), c("yes", "no")))
# 2.25 is the odds ratio of A against B for "yes"
relchk(tab, 2.25, ref = "B", event = "yes")
#> relchk: reported 2.25, to 2 decimals
#> claim: or, ref B, event yes
#> verdict: the claim holds
#>   grp type group ref event value claim
#>  rows   or     B   A    no  2.25 FALSE
#>  rows   or     A   B   yes  2.25  TRUE
#>  cols   or    no yes     B  2.25 FALSE
#>  cols   or   yes  no     A  2.25 FALSE
# 0.44 is the same contrast read with the reference group swapped
relchk(tab, 0.44, ref = "B", event = "yes")
#> relchk: reported 0.44, to 2 decimals
#> claim: or, ref B, event yes
#> verdict: only readings the claim does not allow match
#>   grp type group ref event     value claim
#>  rows   or     B   A   yes 0.4444444 FALSE
#>  rows   or     A   B    no 0.4444444 FALSE
#>  cols   or    no yes     A 0.4444444 FALSE
#>  cols   or   yes  no     B 0.4444444 FALSE
# a risk ratio reported as an odds ratio
relchk(tab, 2.0, type = "or", digits = 1)
#> relchk: reported 2, to 1 decimal
#> claim: or
#> verdict: only readings the claim does not allow match
#>   grp type group ref event value claim
#>  rows   rr     A   B   yes     2 FALSE
```
