# An event study: the effect period by period around adoption

Leads and lags of a policy's adoption from a two-way fixed-effects
regression on indicators of the time since adoption, with the joint test
that the leads are zero (no effect before the policy, the testable part
of parallel trends), and the raw mean outcome of each adoption cohort
around its adoption.

## Usage

``` r
evtst(
  y,
  unit,
  time,
  first,
  leads = 4,
  lags = 4,
  ref = -1,
  bin = TRUE,
  adj = NULL,
  clus = NULL,
  w = NULL,
  ma = 1,
  norm = TRUE,
  mbar = NULL,
  target = NULL,
  draws = 5000,
  seed = NULL,
  level = 0.95
)
```

## Arguments

- y:

  The outcome: finite numbers, one per row (counts from 0 for
  `fam = "pois"` or `"nb"`).

- unit:

  The unit of each row: labels.

- time:

  The period of each row: whole numbers (a year, a month counted from a
  start); each unit at most once per period.

- first:

  The period each row's unit is first treated: numbers, constant within
  units; `NA`, 0 or `Inf` for a unit never treated. Treatment is taken
  to stay on.

- leads, lags:

  The periods before and after adoption estimated: whole numbers,
  `leads` at least 1 and `lags` at least 0.

- ref:

  The event time left out: a whole number from `-leads` to `lags`.

- bin:

  `TRUE` to bin event times beyond the window into its ends, `FALSE` to
  drop their rows.

- adj:

  `NULL` or time-varying controls: a data frame or list of columns, a
  numeric matrix, or one vector, one value per row. Factor, text and
  logical columns enter as indicators of each level but the first.

- clus:

  `NULL` (cluster by unit), the cluster of each row (labels), or a list
  or data frame of two such columns for two-way clustering.

- w:

  `NULL` or weights from 0, one per row (population, matching weights);
  rows of weight 0 take no part.

- ma:

  The width of the moving average of the cohort means: a whole number,
  at least 1 (1 for none).

- norm:

  `TRUE` to scale each cohort's means to 1 in the period before
  adoption.

- mbar:

  `NULL`, or the bounds on the departure from parallel trends relative
  to the largest before adoption: numbers from 0.

- target:

  `NULL` for the mean effect after adoption, or one event time from 0.

- draws:

  The normal draws of the sensitivity interval: a whole number from 100.

- seed:

  The seed of the draws: `NULL` to draw it from R's generator, or a
  whole number from 0.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `evtst`, a list: `est`, one row per event time
(`evtm`, `est`, `se`, `lo`, `hi`, `stat`, `p`), the reference a row with
`est` 0; `pre`, the pre-trend test (`f`, `df1`, `df2`, `p`; `NA` without
leads to test); `vcov`, the variance of the estimated event times;
`coh`, the cohort means (`cohort`, `evtm`, `value`); `n`, the rows used;
`ncl`, the clusters; and `level`; with `mbar`, `sens` (`mbar`, `lob`,
`upb`, `lo`, `hi`) and `brk`, else `NULL`.

## Details

With `G_i` the period unit `i` is first treated and `k = t - G_i`, the
model is
`y_it = sum_k g_k 1(t - G_i = k) + adj_it'th + a_i + m_t + e_it` for `k`
from `-leads` to `lags` but `ref`, the period left out (the one before
adoption by default). Times beyond the window are binned into its ends
(`bin = TRUE`) or their rows dropped. Units never treated have every
indicator 0. The variance, the clusters and the degrees of freedom are
those of
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md).
The pre-trend test is the Wald statistic `W = g_pre' V_pre^-1 g_pre`
over the leads, reported as `F = W / q` on `q` and `G - 1` degrees of
freedom.

With staggered adoption and effects that change over time, the lead and
lag coefficients can mix in effects from other periods (Sun and Abraham
2021);
[`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md)
gives event-time effects free of that.

`coh` is the mean outcome of each cohort (the units first treated in one
period) in each period, as a trailing moving average over `ma` periods,
by time since adoption within the window, and with `norm = TRUE` scaled
to 1 in the period before adoption: a check, before any model, that the
cohorts moved together before the policy.

With `mbar`, the effects are bounded when the trend after adoption may
depart from parallel by up to `mbar` times the largest departure seen
before it (the relative magnitudes of Rambachan and Roth 2023): with `D`
the largest change between consecutive pre-period coefficients (the
reference at 0), the bias at event time `e >= 0` is at most
`(e + 1) mbar D`, so the identified set of the target (the mean effect
after adoption, or the event time `target`) is its estimate `+-` the
mean of `(e + 1)` times `mbar D` (`lob`, `upb`). The interval (`lo`,
`hi`) draws the coefficients `draws` times from their normal
distribution with the cluster-robust covariance, computes the bounds for
each, and takes the `(1 - level) / 2` quantile of the lower and the
`(1 + level) / 2` quantile of the upper: by Bonferroni it covers the
identified set with probability at least `level`. It is conservative;
the conditional and hybrid intervals of Rambachan and Roth are shorter.
`brk` is the smallest `mbar` (by bisection, up to 10; `Inf` beyond)
whose interval holds 0.

## References

Sun, L. and Abraham, S. (2021). Estimating dynamic treatment effects in
event studies with heterogeneous treatment effects. *Journal of
Econometrics*, 225(2), 175-199.
[doi:10.1016/j.jeconom.2020.09.006](https://doi.org/10.1016/j.jeconom.2020.09.006)

Freyaldenhoven, S., Hansen, C. and Shapiro, J. M. (2019). Pre-event
trends in the panel event-study design. *American Economic Review*,
109(9), 3307-3338.
[doi:10.1257/aer.20180609](https://doi.org/10.1257/aer.20180609)

## Examples

``` r
d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
i <- match(d$town, LETTERS)
d$first <- c(2019, 2019, 2019, 2019, 2021, 2021, 2021, 2021,
             NA, NA, NA, NA)[i]
on <- !is.na(d$first) & d$year >= d$first
d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * on + sin(seq_len(120))
e <- evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3)
e
#> evtst: 120 rows, 12 clusters
#>  evtm         est        se         lo         hi        stat            p
#>    -3 -0.09567324 0.5231309 -1.2470766  1.0557301  -0.1828859 8.582148e-01
#>    -2  0.05394300 0.3573220 -0.7325173  0.8404033   0.1509647 8.827364e-01
#>    -1  0.00000000        NA         NA         NA          NA           NA
#>     0 -4.13423655 0.2277190 -4.6354427 -3.6330304 -18.1549928 1.505547e-09
#>     1 -4.25980824 0.4295541 -5.2052505 -3.3143660  -9.9168141 8.036081e-07
#>     2 -4.29625671 0.6385691 -5.7017379 -2.8907756  -6.7279431 3.251955e-05
#>     3 -4.13228464 0.4571512 -5.1384677 -3.1261016  -9.0392074 2.008834e-06
#> leads jointly zero: F = 0.07293073 on 2 and 11 df, p = 0.9301109
e$pre
#>            f df1 df2         p
#> 1 0.07293073   2  11 0.9301109
# how far trends may bend before the effect is lost
evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3,
      mbar = c(0, 1, 2), seed = 1)$sens
#>   mbar       lob       upb        lo         hi
#> 1    0 -4.205647 -4.205647 -4.930333 -3.4716191
#> 2    1 -4.579687 -3.831606 -7.339214 -1.5703317
#> 3    2 -4.953728 -3.457565 -9.831768  0.8305316
```
