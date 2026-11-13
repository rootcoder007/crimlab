# Spatial regression: lag, error and combined models

Fits the spatial lag model `y = rho W y + X beta + e`, the spatial error
model `y = X beta + u` with `u = lambda W u + e`, or both at once
(`"sac"`, one weights matrix for both), with normal errors, by maximum
likelihood (Ord 1975; Anselin 1988): the likelihood is concentrated in
`rho` or `lambda` and maximised by Brent's method on `(-1 / m, 1 / m)`,
`m` the largest row sum of `W` (which bounds its spectral radius), with
\`log \|I

- rho
  W\|`from a dense LU factorisation, so the cost grows as`n^3`per evaluation (at most 5000 units). The combined model is maximised over both parameters by Nelder and Mead's simplex. Standard errors come from the inverse of the expected information matrix (Anselin 1988, ch. 6); these are the estimates of`spatialreg::lagsarlm()`, `errorsarlm()`and`sacsarlm()`. `durbin
  = TRUE`adds the lags`W
  x`of the varying columns (the spatial Durbin model, or the Durbin error model).`method
  =
  "2sls"`fits the lag model by two-stage least squares with the instruments`W
  x`and`W^2
  x`(Kelejian and Prucha 1998), as`spatialreg::stsls()`does, without a likelihood; its standard errors divide the residual sum of squares by`n\`
  less the coefficients.

## Usage

``` r
sarml(y, x, w, model = "lag", durbin = FALSE, method = "ml", int = TRUE)
```

## Arguments

- y:

  The outcome: finite numbers, one per unit.

- x:

  The covariates: a numeric vector, matrix or data frame with one row
  per unit (no intercept column when `int = TRUE`).

- w:

  The weights:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` object, or an `n x n` matrix of weights from 0
  (zero diagonal).

- model:

  `"lag"`, `"error"` or `"sac"`.

- durbin:

  `TRUE` to add the lagged covariates.

- method:

  `"ml"` or `"2sls"` (the lag model only).

- int:

  `TRUE` to add an intercept.

## Value

A list of class `sarml`: `coef`, a data frame of `term` (the covariates,
`lag.` before a lagged one, then `rho` and/or `lambda`), `est`, `se`,
`z`, `p` (two-sided normal), `ols` and `olsse`; `fit`, a named vector:
`loglik` and `aic` (`NaN` for two-stage least squares), `s2` (the error
variance, the residual sum of squares over `n`), `llols` (the
log-likelihood of ordinary least squares) and `n`; `vcov`, the
covariance matrix of the coefficients, rows and columns for every
covariate, `rho` and `lambda` (`NaN` for the one not in the model);
`model`, `method`, `durbin`, `nlag` (the lagged covariates) and `w`.

## Details

Next to each coefficient is ordinary least squares on the same columns
and its usual standard error: leaving out a real spatial lag biases
those coefficients, leaving out a spatial error leaves them unbiased but
their standard errors wrong. In a lag model a coefficient is not the
effect of its variable: see
[`spill()`](https://rootcoder007.github.io/crimlab/reference/spill.md).

## References

Anselin, L. (1988). *Spatial Econometrics: Methods and Models*. Kluwer,
Dordrecht.
[doi:10.1007/978-94-015-7799-1](https://doi.org/10.1007/978-94-015-7799-1)

Kelejian, H. H. and Prucha, I. R. (1998). A generalized spatial
two-stage least squares procedure for estimating a spatial
autoregressive model with autoregressive disturbances. *Journal of Real
Estate Finance and Economics*, 17(1), 99-121.
[doi:10.1023/A:1007707430416](https://doi.org/10.1023/A%3A1007707430416)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(80), y = runif(80))
w <- spwt(p, k = 5)
W <- matrix(0, 80, 80)
W[cbind(w$i, w$j)] <- w$w
x <- rnorm(80)
y <- solve(diag(80) - 0.5 * W, 1 + 2 * x + rnorm(80))
f <- sarml(y, x, w)
f
#> sarml: lag model by ml, 80 units, log-likelihood -119.34 (OLS -128.625)
#>         term       est         se         z            p      ols     olsse
#>  (Intercept) 1.0772273 0.20235921  5.323342 1.018778e-07 1.973029 0.1375247
#>            x 1.9786039 0.12851286 15.396156 1.736932e-53 1.900068 0.1489361
#>          rho 0.4874943 0.08976277  5.430919 5.606464e-08      NaN       NaN
spill(f)
#>   term   direct indirect    total       sed       sei       set
#> 1    x 2.084212 1.776436 3.860648 0.1485669 0.6471067 0.7425332
```
