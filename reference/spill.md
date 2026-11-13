# Direct, indirect and total effects of a spatial lag model

In a model with a spatial lag of the outcome,
`y = rho W y + X beta + W X theta + e`, a change in covariate `k` at
every unit moves the outcomes by
`S_k = (I - rho W)^-1 (beta_k I + theta_k W)`: the mean of its diagonal
is the direct effect, the mean row sum the total, and their difference
the indirect effect, the spillover to other units (LeSage and Pace
2009). An analysis that reads `beta_k` alone can report "no effect"
where the effect is all spillover. With row-standardised weights the
total is `(beta_k + theta_k) / (1 - rho)`, proved in Lean: see
`proof("spill")`. Standard errors are by the delta method in `rho`,
`beta_k` and `theta_k`, with the exact derivatives
`d (I - rho W)^-1 / d rho = (I - rho W)^-1 W (I - rho W)^-1`. The
effects are those of `spatialreg::impacts()`, which gives simulated
rather than delta-method standard errors.

## Usage

``` r
spill(fit = NULL, w = NULL, rho = NULL, b = NULL, th = NULL)
```

## Arguments

- fit:

  A
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  fit of a lag or combined model (`"lag"` or `"sac"`, with or without
  `durbin`); or `NULL`, with `w`, `rho`, `b` (and `th`) given.

- w:

  The weights:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` object, or an `n x n` matrix of weights from 0
  (zero diagonal).

- rho:

  The spatial lag coefficient.

- b:

  The coefficients of the covariates (named, or `x1`, `x2`, ...).

- th:

  `NULL`, or the coefficients of their lags, as many as `b`.

## Value

A data frame with one row per covariate: `term`, `direct`, `indirect`,
`total` and their standard errors `sed`, `sei`, `set` (`NaN` without a
fit).

## References

LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial
Econometrics*. Chapman and Hall/CRC, Boca Raton.
[doi:10.1201/9781420064254](https://doi.org/10.1201/9781420064254)

## Examples

``` r
set.seed(2)
p <- cbind(x = runif(40), y = runif(40))
w <- spwt(p, k = 4)
spill(w = w, rho = 0.4, b = c(gent = 0.1), th = 0.3)
#>   term    direct  indirect     total sed sei set
#> 1 gent 0.1373729 0.5292938 0.6666667 NaN NaN NaN
```
