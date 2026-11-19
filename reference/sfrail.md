# Spatial frailty and cure models

`sfrail()` fits the Weibull proportional hazards model with
log-frailties of the areas under an intrinsic conditional
autoregression, `h(t) = rho t^(rho - 1) exp(x` `beta + u_area)`
(Banerjee, Wall and Carlin 2003); `sfcure()` the promotion-time cure
model (Chen, Ibrahim and Sinha 1999), `S(t) = exp(-theta (1 -`
`exp(-G(t))))`, `theta = exp(x beta + u_area)` and
`G(t) = exp(c + rho log t)`, whose cure fraction is `exp(-theta)`. Both
maximise the penalised likelihood, the penalty `tau u'Qu / 2` with `Q`
the Laplacian of the neighbour weights, by Newton-Raphson from 0 on the
analytic gradient with a Hessian by central differences, the step halved
while the objective falls, until a step is below `tol` (200 iterations
at most). The standard errors come from the last Hessian; the first is
that of the log shape.

## Usage

``` r
sfrail(time, event, x = NULL, reg, w, tau = 1, tol = 1e-10)

sfcure(time, event, x = NULL, reg, w, tau = 1, tol = 1e-10)
```

## Arguments

- time:

  The times: positive numbers.

- event:

  The events: 1 (or `TRUE`) for a failure, 0 for a censored time.

- x:

  `NULL`, or the covariates (no intercept: the area effects hold the
  level), a numeric matrix (or vector) with one row per subject.

- reg:

  The area of each subject: whole numbers from 1 to the number of areas
  of `w`.

- w:

  The neighbours:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` or an `n x n` matrix of weights from 0 (its
  non-zero entries off the diagonal).

- tau:

  The precision of the area effects: one number from 0.

- tol:

  The convergence tolerance: one number strictly between 0 and 1.

## Value

For `sfrail()`, a list of class `sfrail`: `pll` (the penalised
log-likelihood), `shape`, `coef`, `frail`, `se`, `iter`, `tau`. For
`sfcure()`, of class `sfcure`: `pll`, `shape`, `lscale`, `coef`,
`effect`, `se`, `cure` (each subject's cure fraction), `iter`, `tau`.
[`summary()`](https://rdrr.io/r/base/summary.html) gives the table of
estimates.

## References

Banerjee, S., Wall, M. M. and Carlin, B. P. (2003). Frailty modeling for
spatially correlated survival data, with application to infant mortality
in Minnesota. *Biostatistics*, 4(1), 123-142.
[doi:10.1093/biostatistics/4.1.123](https://doi.org/10.1093/biostatistics/4.1.123)

Chen, M.-H., Ibrahim, J. G. and Sinha, D. (1999). A new Bayesian model
for survival data with a surviving fraction. *Journal of the American
Statistical Association*, 94(447), 909-919.
[doi:10.1080/01621459.1999.10474196](https://doi.org/10.1080/01621459.1999.10474196)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
# simulated survival times in three areas
w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
i <- 0:29
x <- sin(i)
t <- 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * x - 0.3 * (i %% 3))
f <- sfrail(t, as.numeric(i %% 5 != 0), x, i %% 3 + 1, w)
f
#> sfrail: Weibull frailty model, shape 3.065, penalised log-likelihood -20.2083
#> [1]  0.9000010 -1.9361730 -1.2253279 -0.8237804
summary(f)
#>       term        est        se         z
#> 1 logshape  1.1201896 0.1648807  6.793941
#> 2       x1  0.9000010 0.3560158  2.527981
#> 3   frail1 -1.9361730 0.4750970 -4.075321
#> 4   frail2 -1.2253279 0.4050962 -3.024783
#> 5   frail3 -0.8237804 0.3803217 -2.166009
```
