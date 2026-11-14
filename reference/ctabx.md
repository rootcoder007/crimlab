# Cross-tabulation with tests and measures of association

Builds a two-way table from labels (counting rows, summing a count
column of an aggregate table, or counting distinct ids per cell) and
gives Pearson's chi-square test, the likelihood-ratio test, Cramer's V
with an interval and its bias-corrected form, the contingency
coefficient, phi for a 2 x 2 table, the share on the diagonal (the rows
whose label is also the column's, as a region at placement against the
region now), and on request Fisher's exact test with, for a 2 x 2 table,
the conditional odds ratio and its exact interval.

## Usage

``` r
ctabx(
  row,
  col = NULL,
  w = NULL,
  agg = NULL,
  correct = FALSE,
  exact = FALSE,
  level = 0.95
)
```

## Arguments

- row:

  The row label of each record (a factor by its levels, text or numbers
  by their sorted values); or a whole table: a numeric matrix or 2-d
  table of counts, with `col` `NULL`; or a `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md).

- col:

  The column label of each record.

- w:

  `NULL`, or the value of each record that `agg` combines into its cell:
  a count (`"sum"`), a measure (`"max"`, `"min"`, `"mean"`), or an id
  (`"dist"`, the distinct ids per cell).

- agg:

  How the records of a cell combine: `"n"` (count them, the default
  without `w`), `"sum"` (the default with `w`), `"max"`, `"min"`,
  `"mean"` or `"dist"`. An empty cell is 0.

- correct:

  `TRUE` for the continuity correction of a 2 x 2 chi-square, as
  [`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)'s
  default; `FALSE` (here the default) for Pearson's statistic.

- exact:

  `TRUE` for Fisher's exact test.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `ctabx`, a list: `tab`, the table (all labels);
`expd` and `res`, the expected counts and Pearson residuals of the table
without empty rows and columns; `stat`, a one-row data frame: `n`,
`chisq`, `df`, `p`, `g2`, `g2p` (the likelihood-ratio statistic and its
p-value), `v`, `vlo`, `vhi`, `vbc` (bias-corrected V), `cc` (the
contingency coefficient), `phi`, `diag`, `fp` (Fisher's p-value), `or`,
`orlo`, `orhi`; and `drop`, `agg`, `correct` and `level`.

## Details

Rows or columns with no counts are dropped before the statistics (their
number is `drop`), as an unused factor level carries no information; the
table `tab` keeps them. Cramer's V is `sqrt(chisq / (n (k - 1)))`, `k`
the smaller dimension, from the uncorrected statistic even when
`correct` is `TRUE`; its interval inverts the noncentral chi-square
distribution, `sqrt(ncp / (n (k - 1)))`, as the default of
`DescTools::CramerV()`; the bias-corrected V is Bergsma's. The exact
test enumerates every table with the observed margins and refuses one
that would take more than 1e7 tables. With `agg` other than `"n"`,
`"sum"` or `"dist"` the table holds values, not counts, and the
statistics are `NA`.

## References

Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley. ISBN
978-0-470-46363-5.

Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
*Journal of the Korean Statistical Society*, 42(3), 323-328.
[doi:10.1016/j.jkss.2012.10.002](https://doi.org/10.1016/j.jkss.2012.10.002)

## Examples

``` r
# force type by weapon, one row per incident
wp <- c("none", "none", "knife", "gun", "none", "knife", "gun", "none")
fc <- c("hands", "hands", "CEW", "firearm", "CEW", "hands", "firearm",
        "hands")
ctabx(wp, fc)
#> ctabx: 3 x 3 table (n), n = 8
#>       CEW firearm hands
#> gun     0       2     0
#> knife   1       0     1
#> none    1       0     3
#> chi-square 8.5 on 4 df, p = 0.07488723
#> Cramer's V 0.728869 [0, 1], bias-corrected 0.586302
# region at placement against region now: the diagonal share
ctabx(c("N", "N", "E", "W", "E"), c("N", "E", "E", "W", "E"))
#> ctabx: 3 x 3 table (n), n = 5
#>   E N W
#> E 2 0 0
#> N 1 1 0
#> W 0 0 1
#> chi-square 6.666667 on 4 df, p = 0.1545873
#> Cramer's V 0.8164966 [0, 1], bias-corrected 0.5773503
#> diagonal share 0.8
# a published aggregate table: sum its count column
d <- data.frame(race = c("A", "A", "B", "B"), sex = c("F", "M", "F", "M"),
                n = c(12, 40, 9, 51))
ctabx(d$race, d$sex, w = d$n, exact = TRUE)
#> ctabx: 2 x 2 table (sum), n = 112
#>    F  M
#> A 12 40
#> B  9 51
#> chi-square 1.192899 on 1 df, p = 0.2747455
#> Cramer's V 0.1032031 [0, 0.2883748], bias-corrected 0.04070387
#> Fisher's exact p = 0.3348875
```
