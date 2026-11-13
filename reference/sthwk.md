# Space-time Hawkes process

Fits the self-exciting point process of crime in space and time of
Mohler et al. (2011): \$\$\lambda(t, x, y) = \mu + \eta \sum\_{t_j \< t}
\beta e^{-\beta (t - t_j)} \phi\_\sigma(x - x_j, y - y_j),\$\$ a
constant background \\\mu\\ per unit area and time over the window and,
after each event, `eta` expected offspring at exponential lags (rate
`beta`) and normal offsets (standard deviation `sigma`). The
log-likelihood takes the spatial kernel's mass as 1 inside the window
(events near its edge lose a little offspring mass outside); it is
maximised with its analytic gradient by BFGS over the logarithms of
`mu`, `beta`, `sigma` and the logit of `eta`, from two starts or
`start`. With `par` nothing is fitted: the log-likelihood is evaluated
there. `at` gives the conditional intensity at query points, from the
events strictly before each.
[`hwksim()`](https://rootcoder007.github.io/crimlab/reference/hwksim.md)
simulates from a fit.

## Usage

``` r
sthwk(
  times,
  x,
  y,
  hor,
  win = NULL,
  par = NULL,
  at = NULL,
  start = NULL,
  maxit = 500
)
```

## Arguments

- times:

  Event times: numbers sorted in `[0, hor]`, at least 10.

- x, y:

  The places of the events: finite numbers in planar units.

- hor:

  The end of the window in time: one positive number.

- win:

  `NULL` (the events' bounding box) or the rectangle
  `c(xmin, xmax, ymin, ymax)` holding every event.

- par:

  `NULL` to fit, or the parameters `mu`, `eta`, `beta`, `sigma` at which
  to evaluate.

- at:

  `NULL` or query points: a data frame or list with `t`, `x` and `y`.

- start:

  `NULL` or the starting parameters, as `par`.

- maxit:

  The most BFGS steps per start: one whole number from 1.

## Value

An object of class `sthwk`, a list: `par`, `se` and `cov` (`mu`, `eta`,
`beta`, `sigma`); `loglik`, `aic`; `conv`, `iter`; `fit` (`FALSE` when
`par` was given); `n`, `hor`, `win`, `area`; and `lam`, the intensities
at `at` (`NULL` without it).

## References

Mohler, G. O., Short, M. B., Brantingham, P. J., Schoenberg, F. P. and
Tita, G. E. (2011). Self-exciting point process modeling of crime.
*Journal of the American Statistical Association*, 106(493), 100-108.
[doi:10.1198/jasa.2011.ap09546](https://doi.org/10.1198/jasa.2011.ap09546)

Reinhart, A. (2018). A review of self-exciting spatio-temporal point
processes and their applications. *Statistical Science*, 33(3), 299-318.
[doi:10.1214/17-STS629](https://doi.org/10.1214/17-STS629)

## Examples

``` r
s <- hwksim(c(log(0.02), 0.5, 1), 100, win = c(0, 10, 0, 10),
            sigma = 0.3, seed = 3)
f <- sthwk(s$t, s$x, s$y, 100, win = c(0, 10, 0, 10))
f
#> sthwk: 343 events in [0, 100] over an area of 100
#>              est          se
#> mu    0.01908445 0.001555915
#> eta   0.45225241 0.042523094
#> beta  0.95205502 0.114791215
#> sigma 0.29276704 0.021369007
#> log-likelihood -1288.096
sthwk(s$t, s$x, s$y, 100, win = c(0, 10, 0, 10), par = f$par,
      at = data.frame(t = 50, x = 5, y = 5))$lam
#> [1] 0.01938062
```
