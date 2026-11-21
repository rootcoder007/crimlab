# Cramer's V of a two-way table

The strength of association in a contingency table as one number: 0 when
the rows and columns are independent in the sample, 1 when each row
holds a single column. It is the measure
[`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
reports beside its tests, given here on its own for a table that is
already counted.

## Usage

``` r
cramerv(tbl, biascor = TRUE, level = 0.95)
```

## Arguments

- tbl:

  A two-way table of counts: a numeric matrix or a 2-d table.

- biascor:

  `TRUE` (the default) for Bergsma's bias-corrected V in `v`; `FALSE`
  for the uncorrected one.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

A named numeric vector of class `cramerv`: `v`, `vraw`, `lo`, `hi`,
`chisq`, `df`, `p`, `n`, `minexp` (the smallest expected count) and
`drop` (the empty rows and columns left out). In a table with fewer than
two rows or columns that hold counts the statistics are `NA`.

## Details

Cramer's V is `sqrt(chisq / (n (k - 1)))`, `k` the smaller dimension,
from Pearson's uncorrected statistic. In a small table it is biased
upward (it cannot be negative, so sampling noise alone gives a positive
value); with `biascor = TRUE`, the default, `v` is Bergsma's
bias-corrected form, and `vraw` holds the uncorrected one. The interval
`lo` to `hi` inverts the noncentral chi-square distribution and belongs
to the uncorrected V. Rows and columns with no counts are dropped first,
as an unused factor level carries no information.

The p-value is that of the chi-square approximation, which rests on the
expected counts: `minexp` is the smallest of them, and below about 5 the
exact test of
[`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
(`exact = TRUE`) is the one to report.

Errors: `tbl` not a numeric matrix or two-way table of finite,
non-negative counts; `biascor` not `TRUE` or `FALSE`; a `level` outside
(0, 1).

## References

Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton
University Press.

Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
*Journal of the Korean Statistical Society*, 42(3), 323-328.
[doi:10.1016/j.jkss.2012.10.002](https://doi.org/10.1016/j.jkss.2012.10.002)

## See also

[`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
for the table from records, the tests and the exact test.

## Examples

``` r
# use of force by weapon, counted
tab <- matrix(c(30, 6, 2, 9, 14, 3, 1, 4, 11), 3,
              dimnames = list(c("none", "knife", "gun"),
                              c("hands", "CEW", "firearm")))
cramerv(tab)
#> cramerv: V = 0.5150056 (bias-corrected; uncorrected 0.5327594, 95% interval 0.3563396 to 0.6720259)
#> chi-square 45.41321 on 4 df, p = 3.262244e-09, n = 80
#> smallest expected count 3.2: report the exact test, ctabx(exact = TRUE)
cramerv(tab, biascor = FALSE)[["v"]]
#> [1] 0.5327594
# independence: the corrected V is 0
cramerv(matrix(c(10, 20, 10, 20), 2))[["v"]]
#> [1] 0
try(cramerv(1:4))
#> Error : cramerv: `tbl` must be a two-way table of counts (a numeric matrix or a 2-d table).
```
