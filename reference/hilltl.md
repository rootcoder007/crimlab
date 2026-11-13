# Tail index of a heavy-tailed count

The exponent `alpha` of a power-law tail, `P(X = k)` proportional to
`k^-alpha` above a threshold `xmin`, as when a few people account for
most repeat placements or most steps of an offender's path are short and
a few very long (Clauset, Shalizi and Newman 2009). For counts the
maximum likelihood estimate of the discrete law is the default, its
normalising constant the Hurwitz zeta function `zeta(alpha, xmin)`;
`approx = TRUE` gives the closed form with `xmin - 1/2`, and
`discrete = FALSE` the continuous Hill estimator
`1 + n / sum log(x / xmin)` for measurements such as distances. `ks` is
the Kolmogorov-Smirnov distance between the tail and the fitted law: an
exponent is only as good as that fit. The values at or below 0 are left
out.

## Usage

``` r
hilltl(x, xmin = NULL, discrete = TRUE, approx = FALSE, mintl = 3)
```

## Arguments

- x:

  The values: numbers from 0.

- xmin:

  `NULL` (the 50th largest value, or the smallest when there are fewer
  than 50) or the threshold: one positive number.

- discrete:

  `TRUE` for counts, `FALSE` for measurements.

- approx:

  `TRUE` for the closed form for counts instead of the maximum
  likelihood estimate.

- mintl:

  The fewest tail values to fit to: a whole number (at least 2 are
  always needed).

## Value

An object of class `hilltl`, a list: `alpha` and its standard error `se`
(the observed information for the discrete estimate,
`(alpha - 1) / sqrt(n)` otherwise; both `NA` when no fit is made),
`xmin`, `ntl` (the tail size), `ks`, `reliable` (`TRUE` with 50 or more
tail values and an estimate inside the search interval `[1.0001, 25]`)
and `method`.

## References

Clauset, A., Shalizi, C. R. and Newman, M. E. J. (2009). Power-law
distributions in empirical data. *SIAM Review*, 51(4), 661-703.
[doi:10.1137/070710111](https://doi.org/10.1137/070710111)

Hill, B. M. (1975). A simple general approach to inference about the
tail of a distribution. *The Annals of Statistics*, 3(5), 1163-1174.
[doi:10.1214/aos/1176343247](https://doi.org/10.1214/aos/1176343247)

## Examples

``` r
set.seed(4)
k <- floor(runif(400)^(-1 / 1.5))
hilltl(k, xmin = 1)
#> hilltl: discrete maximum likelihood
#> alpha 2.146 (se 0.0619), xmin 1, 400 tail values, KS 0.0377
```
