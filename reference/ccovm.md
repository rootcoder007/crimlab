# Covariance matrix (C backend)

The covariance of every pair of columns of a numeric matrix or data
frame. Column names are carried to both dimensions of the result.

## Usage

``` r
ccovm(x)
```

## Arguments

- x:

  A numeric matrix, or a data frame of numeric columns (rows
  observations, columns variables).

## Value

A symmetric `ncol(x)` by `ncol(x)` matrix.

## Details

Entry \\(j, k)\\ is \\\frac{1}{n - 1} \sum_i (x\_{ij} - \bar
x_j)(x\_{ik} - \bar x_k)\\, as
[`stats::cov()`](https://rdrr.io/r/stats/cor.html) gives it; the
diagonal holds the column variances. Missing values are not removed: an
`NA` in a column makes every covariance with that column `NA`. A single
row gives `NaN` throughout. A plain vector, and a data frame with a
column that is not numeric, are refused.

## See also

[`ccorr()`](https://rootcoder007.github.io/crimlab/reference/cstat.md),
[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md).

## Examples

``` r
X <- cbind(a = c(1, 2, 3, 4), b = c(2, 4, 7, 8), c = c(5, 3, 2, 1))
ccovm(X)
#>           a         b         c
#> a  1.666667  3.500000 -2.166667
#> b  3.500000  7.583333 -4.583333
#> c -2.166667 -4.583333  2.916667
all.equal(ccovm(X), stats::cov(X))
#> [1] TRUE
ccovm(data.frame(u = 1:5, v = c(2, 1, 4, 3, 6)))
#>     u   v
#> u 2.5 2.5
#> v 2.5 3.7
try(ccovm(1:5))
#> Error : ccovm: `x` must be a numeric matrix or data frame.
```
