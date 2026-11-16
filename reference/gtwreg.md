# Geographically and temporally weighted regression

The regression of Kuo and Putra (2021, eqs 4-9), after Huang, Wu and
Barry (2010): at each observation `i` the least-squares fit of all
observations weighted by a Gaussian kernel of their space-time distance,
so that the coefficients vary over places and periods. They relate
robbery, domestic violence and fraud in the local government areas of
New South Wales to the air quality index month by month.

## Usage

``` r
gtwreg(
  y,
  x,
  u,
  v,
  t = NULL,
  h,
  tau = 0,
  crit = "cv",
  lonlat = FALSE,
  cons = TRUE
)
```

## Arguments

- y:

  The outcome: finite numbers.

- x:

  The regressors: a numeric matrix or data frame, one row per value of
  `y` (columns named), or `NULL` for the constant alone.

- u, v:

  The coordinates of each observation.

- t:

  `NULL` or the time of each observation (numbers).

- h:

  The bandwidths tried: positive numbers.

- tau:

  The ratios of the time to the space scale tried: numbers from 0.

- crit:

  `"cv"` or `"aicc"`: how `h` and `tau` are chosen.

- lonlat:

  `TRUE` for longitude and latitude in degrees.

- cons:

  `TRUE` for a local constant.

## Value

An object of class `gtwreg`, a list: `coef` and `se`, matrices of the
local coefficients and their standard errors (one row per observation);
`fit`, `resid`, `hat` and `loo` (the leave-one-out residuals); `stat`, a
named vector (`h`, `tau`, `rss`, `trs`, `trss`, `sigma2`, `cv`, `aicc`,
`edf`, `r2`); and `grid`, the criteria of every pair tried.

## Details

With the spatial distance `ds` (planar, or kilometres on the sphere with
`lonlat`) and the time distance `dt`, the space-time distance is
`d^2 = ds^2 + tau dt^2` (eq 7, the spatial scale set to 1) and the
weight `w = exp(-d^2 / h^2)` (eq 8); `beta_i = (X'W_i X)^-1 X'W_i y` (eq
5). Without `t` (or with `tau = 0`) it is geographically weighted
regression. The cross-validation score is `sum (y_i - y_(-i))^2` (eq 9),
the leave-one-out prediction from the fit at `i` without observation `i`
(computed as `e_i / (1 - S_ii)`, `S` the hat matrix). Kuo and Putra
choose `tau` and the bandwidth `h` by cross-validation or by "the
corrected Akaike information criterion" without giving it; the `aicc`
here is the one of the GWR literature,
`n log(RSS / n) + n log(2 pi) + n (n + tr S) / (n - 2 - tr S)`. Give
several values of `h` and `tau` and the pair with the smallest `crit` is
fitted; all are in `grid`. The local standard errors are `sigma^2` times
the diagonal of `(X'W_iX)^-1 X'W_i^2 X (X'W_iX)^-1`,
`sigma^2 = RSS / (n - 2 tr S + tr S'S)`.

## References

Kuo, P.-F. and Putra, I. G. B. (2021). Analyzing the relationship
between air pollution and various types of crime. *PLoS ONE*, 16(8),
e0255653.
[doi:10.1371/journal.pone.0255653](https://doi.org/10.1371/journal.pone.0255653)

Huang, B., Wu, B. and Barry, M. (2010). Geographically and temporally
weighted regression for modeling spatio-temporal variation in house
prices. *International Journal of Geographical Information Science*,
24(3), 383-401.
[doi:10.1080/13658810802672469](https://doi.org/10.1080/13658810802672469)

## See also

[`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md),
[`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)

## Examples

``` r
set.seed(4)
d <- expand.grid(u = 1:6, v = 1:6, t = 1:4)
d$aqi <- rnorm(nrow(d))
d$y <- 1 + (0.2 + 0.1 * d$u) * d$aqi + rnorm(nrow(d), 0, 0.3)
r <- gtwreg(d$y, d["aqi"], d$u, d$v, d$t, h = c(1.5, 3, 6),
            tau = c(0, 1))
r
#> gtwreg: geographically and temporally weighted regression, h = 1.5, tau = 0
#> CV 13.488, AICc 70.291, R2 0.796, effective parameters 13.35
#> local coefficients (quartiles):
#>                0%       25%       50%       75%      100%
#> (const) 0.9361268 0.9679731 0.9848207 0.9955117 1.0355406
#> aqi     0.2981053 0.4329041 0.5922378 0.6711381 0.7942628
r$grid
#>     h tau       cv     aicc
#> 1 1.5   0 13.48802 70.29134
#> 2 3.0   0 14.16615 75.47577
#> 3 6.0   0 15.53506 88.76016
#> 4 1.5   1 13.72343 74.63537
#> 5 3.0   1 14.16584 75.13253
#> 6 6.0   1 15.53053 88.64167
```
