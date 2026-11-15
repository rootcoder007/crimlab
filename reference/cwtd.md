# Weighted mean and variance (C backend)

The weighted mean and variance of a numeric vector, for data where some
observations count for more than others: an average over areas weighted
by population, or measurements weighted by their precision.

## Usage

``` r
cwtd(x, w)
```

## Arguments

- x:

  Numeric vector of observations.

- w:

  Non-negative finite weights, the length of `x`, not all zero.

## Value

A named numeric vector: `mean` and `variance`.

## Details

The mean is \\\bar x_w = \sum_i w_i x_i / \sum_i w_i\\. The weights are
read as RELIABILITY weights, so the variance is \$\$s_w^2 = \frac{\sum_i
w_i (x_i - \bar x_w)^2} {\sum_i w_i - \sum_i w_i^2 / \sum_i w_i},\$\$
the unbiased estimator under that reading; with every weight 1 it is
[`stats::var()`](https://rdrr.io/r/stats/cor.html). Frequency weights (a
weight of 3 meaning three identical rows) call for \\\sum w - 1\\
instead, and survey design weights for a design-based variance: do not
use this variance for those. An `NA` in either vector gives `NaN` for
both entries.

## See also

[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md),
[`cipw()`](https://rootcoder007.github.io/crimlab/reference/cipw.md).

## Examples

``` r
x <- c(10, 20, 30, 40)
w <- c(1, 1, 2, 4)
cwtd(x, w)
#>     mean variance 
#>  31.2500 169.0476 
all.equal(cwtd(x, w)[["mean"]], stats::weighted.mean(x, w))
#> [1] TRUE
all.equal(cwtd(x, rep(1, 4))[["variance"]], stats::var(x))
#> [1] TRUE
try(cwtd(x, c(1, -1, 1, 1)))
#> Error : cwtd: `w` must be non-negative and finite, not all zero.
```
