# Standard deviation and Euclidean distance (C backend)

`csdev()` returns the standard deviation of a numeric vector and
`cdistm()` the straight-line (Euclidean) distance between two vectors of
the same length.

## Usage

``` r
csdev(x, ddof = 1L)

cdistm(a, b)
```

## Arguments

- x, a, b:

  Numeric (or logical) vectors; `a` and `b` of the same length.

- ddof:

  The delta degrees of freedom, one whole number from 0: the divisor is
  `n - ddof`.

## Value

A number: the standard deviation in the units of `x`, or the distance in
the units of `a` and `b`.

## Details

`csdev()` computes \$\$s = \sqrt{\frac{1}{n - \mathrm{ddof}}
\sum\_{i=1}^n (x_i - \bar x)^2},\$\$ so `ddof = 1` is the sample
standard deviation of [`stats::sd()`](https://rdrr.io/r/stats/sd.html)
and `ddof = 0` the population figure. `cdistm()` computes \\\sqrt{\sum_i
(a_i - b_i)^2}\\.

Missing values propagate: one `NA` makes the result `NA`. A vector with
no more than `ddof` values has no defined spread and gives `NaN`. Two
empty vectors are at distance 0. Text, a factor, a list or a data frame
is refused, and so are `a` and `b` of different lengths; logical input
counts as 0 and 1.

## See also

[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md),
[`cmad()`](https://rootcoder007.github.io/crimlab/reference/crob.md).

## Examples

``` r
csdev(c(2, 4, 4, 4, 5, 5, 7, 9))
#> [1] 2.13809
all.equal(csdev(1:10), stats::sd(1:10))
#> [1] TRUE
csdev(1:10, ddof = 0)
#> [1] 2.872281
csdev(c(3, NA, 5))
#> [1] NA
cdistm(c(0, 0), c(3, 4))
#> [1] 5
try(cdistm(1:3, 1:2))
#> Error : cdistm: `a` and `b` must have the same length.
```
