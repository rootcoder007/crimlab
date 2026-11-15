# Regularized incomplete gamma function (C backend)

The lower regularized incomplete gamma function \\P(a, x)\\, the
distribution function of a gamma distribution with shape `a` and unit
rate; the p-value of a chi-square statistic \\q\\ on \\k\\ degrees of
freedom is \\1 - P(k/2, q/2)\\.

## Usage

``` r
cgamcdf(shape, x)
```

## Arguments

- shape:

  The shape \\a\\, one positive number.

- x:

  Numeric vector of quantiles.

## Value

A numeric vector of probabilities, the length of `x`.

## Details

\$\$P(a, x) = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} \\ dt,\$\$
computed by the package's own gamma kernel, which agrees with
`stats::pgamma(x, shape)` to rounding. For a rate \\r\\, pass `x * r`. A
negative `x` gives 0, `NA` gives `NA` and `Inf` gives 1; `shape = Inf`
gives 0 for every finite `x`, as
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html) does.

## See also

[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html).

## Examples

``` r
cgamcdf(3.5, c(0.5, 1, 4, 12))
#> [1] 0.005171463 0.040159631 0.667406097 0.998860649
all.equal(cgamcdf(3.5, c(0.5, 1, 4, 12)),
          stats::pgamma(c(0.5, 1, 4, 12), shape = 3.5))
#> [1] TRUE
1 - cgamcdf(2 / 2, 5.99 / 2)
#> [1] 0.05003663
try(cgamcdf(0, 1))
#> Error : cgamcdf: `shape` must be one positive number.
```
