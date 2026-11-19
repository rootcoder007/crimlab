# Generalized synthetic control with interactive fixed effects

The effect of a treatment adopted by some units of a balanced panel at
possibly different times, estimated as the generalized synthetic control
of Xu (2017): a model of the untreated outcome fitted to the units never
treated, then carried to each treated unit through its periods before
adoption.

## Usage

``` r
ifesc(
  data,
  y,
  d,
  unit,
  time,
  x = NULL,
  r = NULL,
  rmax = 5,
  se = c("boot", "none"),
  nboot = 200,
  seed = NULL,
  level = 0.95,
  tol = 1e-08,
  maxit = 1000
)

# S3 method for class 'ifesc'
print(x, ...)

# S3 method for class 'ifesc'
summary(object, ...)
```

## Arguments

- data:

  A data frame, one row per unit and time (balanced: every unit at every
  time once).

- y:

  The name of the outcome column (numbers).

- d:

  The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`), staying
  1 once it turns 1; never-treated units are the controls.

- unit, time:

  The names of the unit column (labels) and the time column (numbers or
  dates; ordered by value).

- x:

  `NULL` or the names of covariate columns (numbers).

- r:

  `NULL` (cross-validation) or the number of factors, a whole number
  from 0.

- rmax:

  The largest number of factors cross-validation tries.

- se:

  `"boot"` (the bootstrap of units) or `"none"`.

- nboot:

  The number of bootstrap draws, at least 2.

- seed:

  `NULL` (the generator's current state) or one whole number.

- level:

  The confidence level of the intervals.

- tol:

  The convergence tolerance of `b`.

- maxit:

  The most alternations of Bai's estimator.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html); for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), passed to
  the base-graphics call.

- object:

  An object of class `ifesc`.

## Value

An object of class `ifesc`, a list: `att`, a one-row data frame (`est`,
`se`, `stat`, `p`, `lo`, `hi`); `eff`, one row per event time `e` (`n`
treated units observed then, and the same columns); `r`; `cv`, a data
frame of `r` and `mspe` (`NULL` when `r` was given); `beta`, the
covariates' coefficients; `cf`, one row per treated unit and time
(`unit`, `time`, `e`, `y`, `y0`, `eff`); `factors`, `T x r`; `draws`,
the bootstrap's `att` then `eff` per draw, a row each (`NULL` for
`se = "none"`); `ntr`, `nco`, `nt`; `iter`, `conv`; `nboot`, `nfail`,
`level`. [`summary()`](https://rdrr.io/r/base/summary.html) returns
`eff`.

## Details

The untreated outcome of unit i at time t is
`y_it = x_it'b + mu + alpha_i` `+ xi_t + lambda_i'f_t + e_it`, the
interactive fixed effects model of Bai (2009): unit and time effects
plus `r` common factors `f_t` with unit loadings `lambda_i`, which let
units follow different trends (a difference in differences assumes
`r = 0`). On the never-treated units the estimator of Bai (2009)
alternates the factors, `sqrt(T)` times the first `r` eigenvectors of
`W W'` for the doubly demeaned residuals `W` net of covariates
(normalised `F'F / T = I`), and `b`, by least squares of the doubly
demeaned outcome less the common component on the doubly demeaned
covariates, until `b` changes by less than `tol`. Each treated unit's
`alpha_i` and `lambda_i` come from least squares of its residuals on
`(1, f_t)` over its periods before adoption, which must number more than
`r + 1`; its counterfactual `y0` is the model's value at every period,
and its effect at a period is `y - y0`. `att` averages the effects of
the treated periods; `eff` averages them by periods since adoption
(event time `e`, 0 the first treated period; the effects before adoption
are the fit's in-sample gaps).

`r = NULL` chooses the number of factors from 0 to `rmax` by the
cross-validation of Xu (2017, sec. 3.1, algorithm 1): each treated
unit's pre-adoption periods are left out one at a time, the loadings
fitted on the rest and the left-out outcome predicted; the `r` with the
smallest mean squared prediction error is kept (`cv`). The search stops
at the first `r` some treated unit has too few periods for.

The standard errors are those of a bootstrap of units, the "simple
nonparametric bootstrap" Xu (2017, sec. 3.2, p. 63) finds valid when the
treated units are many: each draw takes as many units as the panel has,
with replacement, and refits with the chosen `r`; draws with no treated
unit or with no more control units than `r` fail and are left out. With
few treated units Xu (2017) prefers a parametric bootstrap that
resamples residuals, which is not built here. Rows of a unit drawn twice
enter twice.

Barber (2026) estimates the effect of the New York City congestion
charge on fine particles, traffic speeds and arrests with the
generalized synthetic control, and the spatial shift of arrests with a
difference in differences (see
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)).

## References

Bai, J. (2009). Panel data models with interactive fixed effects.
*Econometrica*, 77(4), 1229-1279.
[doi:10.3982/ECTA6135](https://doi.org/10.3982/ECTA6135)

Xu, Y. (2017). Generalized synthetic control method: causal inference
with interactive fixed effects models. *Political Analysis*, 25(1),
57-76. [doi:10.1017/pan.2016.2](https://doi.org/10.1017/pan.2016.2)

Barber, T. J. (2026). *The New York City congestion charge: air quality,
traffic, and crime*. MS thesis, Cornell University.

## Examples

``` r
set.seed(1)
p <- expand.grid(id = 1:30, yr = 1:12)
f <- cumsum(rnorm(12))
l <- rnorm(30)
p$tr <- as.integer(p$id <= 6 & p$yr >= 9)
p$y <- l[p$id] * f[p$yr] + 2 * p$tr + rnorm(nrow(p), 0, 0.3)
r <- ifesc(p, "y", "tr", "id", "yr", r = 1, nboot = 50, seed = 1)
r
#> ifesc: generalized synthetic control, 6 treated and 24 control units, 12 times, 1 factor(s)
#>    est     se  stat         p    lo    hi
#>  1.994 0.2222 8.974 2.851e-19 1.558 2.429
```
