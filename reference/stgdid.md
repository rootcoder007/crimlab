# Effects of a policy adopted at different times

Difference-in-differences estimators for staggered adoption that never
use an already-treated unit as a control: the group-time effects of
Callaway and Sant'Anna (2021), optionally with covariates; the
regression on cohort-period cells of Wooldridge (2021); and, for a
treatment that switches on and off, the estimator of de Chaisemartin and
D'Haultfoeuille (2020) with its placebo. Standard errors come from the
cluster bootstrap.

## Usage

``` r
stgdid(
  y,
  unit,
  time,
  first = NULL,
  trt = NULL,
  meth = "csa",
  ctrl = "notyet",
  adj = NULL,
  est = "dr",
  clus = NULL,
  reps = 499,
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

  The period of each row: numbers (a year, a month index); each unit at
  most once per period.

- first:

  The period each row's unit is first treated (constant within units;
  `NA`, 0 or `Inf` for never), or `NULL` with `trt`.

- trt:

  The treatment of each row, 0 or 1, or `NULL` with `first`. For
  `meth = "csa"` and `"wool"` it must stay on once on.

- meth:

  `"csa"`, `"wool"` or `"dcdh"`.

- ctrl:

  `"notyet"` or `"never"`: the comparison units of `"csa"`.

- adj:

  `NULL` or covariates for `"csa"`: a data frame or list of columns, a
  numeric matrix or one vector, one value per row.

- est:

  With `adj`: `"dr"`, `"ipw"` or `"reg"`.

- clus:

  `NULL` (resample units) or the cluster of each row, constant within
  units.

- reps:

  The bootstrap samples: a whole number from 0 (no standard errors).

- seed:

  The seed of the bootstrap: `NULL` to draw it from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it), or a
  whole number from 0.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `stgdid`, a list: `att`, the group-time effects
(`cohort`, `period`, `evtm`, `est`, `se`, `lo`, `hi`, `ncoh` the
cohort's size; `NULL` for `"dcdh"`); `overall` (`term`, `est`, `se`,
`lo`, `hi`, and for `"dcdh"` the rows `"effect"` and `"placebo"` with
`n`, the switchers); `dynamic`, by event time (`evtm`, `est`, `se`,
`lo`, `hi`; `NULL` for `"dcdh"`); `meth`; `ctrl`; `reps`; and `level`.

## Details

Let cohort `g` be the units first treated in period `g`, and the
comparison units those never treated (`ctrl = "never"`) or not yet
treated by period `max(t, b)` (`"notyet"`, outside cohort `g`). For
`meth = "csa"`, `ATT(g, t)` is the cohort's mean change in the outcome
from the base period `b` less the comparison units', with `b = g - 1`
from adoption on and `b = t - 1` before (the `"varying"` base of the
'did' package), so an effect before adoption is a placebo over one
period. With `adj`, parallel trends need hold only given the covariates,
taken at the base period, and each cell is the estimator of Sant'Anna
and Zhao (2020) on the cohort and its comparison units: `"dr"`, the
doubly robust `mean_D(dy - m(x)) - sum w (dy - m(x)) / sum w` with `m`
the least-squares fit of the change among the comparison units and
`w = p / (1 - p)` from the logistic propensity `p` of belonging to the
cohort; `"ipw"` without `m`; `"reg"` without the weights. Comparison
units with `p` of 0.995 or more are left out, as 'DRDID' does. Use
covariates fixed before treatment.

For `meth = "wool"`, the regression of the outcome on unit and period
effects and one indicator for every treated cohort-period cell gives
each `ATT(g, t)` after adoption as a coefficient (without covariates
this is the imputation estimator); without units never treated, the
periods from the last adoption on are left out.

The overall effect averages the `ATT(g, t)` after adoption weighted by
cohort size (the `"simple"` aggregation of 'did'); the effect at event
time `e` averages `ATT(g, g + e)` over the cohorts observed at `e`,
weighted by cohort size. Event time counts observed periods. Units
treated from the first period are left out.

For `meth = "dcdh"`, in each period the units that switch into treatment
are compared with the units untreated in both periods, and the units
that switch out with the units treated in both:
`DID_M = sum_t (N_join,t DID_join,t + N_leave,t DID_leave,t) / N_S` over
the `N_S` switchers that have a comparison. The placebo runs the same
comparison one period earlier, for units whose treatment did not change
over the two periods before; it should be near zero.

Standard errors are the standard deviations over `reps` bootstrap
samples of the clusters (units by default), drawn with replacement, with
everything recomputed; intervals are normal.

## References

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-differences
with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)

Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
difference-in-differences estimators. *Journal of Econometrics*, 219(1),
101-122.
[doi:10.1016/j.jeconom.2020.06.003](https://doi.org/10.1016/j.jeconom.2020.06.003)

Wooldridge, J. M. (2021). Two-way fixed effects, the two-way Mundlak
regression, and difference-in-differences estimators.
[doi:10.2139/ssrn.3906345](https://doi.org/10.2139/ssrn.3906345)

de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed
effects estimators with heterogeneous treatment effects. *American
Economic Review*, 110(9), 2964-2996.
[doi:10.1257/aer.20181169](https://doi.org/10.1257/aer.20181169)

## Examples

``` r
d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
i <- match(d$town, LETTERS)
d$first <- c(2019, 2019, 2019, 2019, 2021, 2021, 2021, 2021,
             NA, NA, NA, NA)[i]
on <- !is.na(d$first) & d$year >= d$first
d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * on + sin(seq_len(120))
s <- stgdid(d$rate, d$town, d$year, d$first, reps = 99, seed = 1)
s$overall
#>   term       est        se        lo        hi
#> 1  att -4.034678 0.4943916 -5.003667 -3.065688
s$dynamic
#>    evtm         est        se         lo         hi
#> 1    -5 -0.37997631 0.3609287 -1.0873837  0.3274310
#> 2    -4 -0.22087536 0.4097379 -1.0239469  0.5821962
#> 3    -3  0.17184724 0.3162810 -0.4480521  0.7917466
#> 4    -2  0.12208141 0.3516011 -0.5670441  0.8112069
#> 5    -1  0.03901677 0.3476043 -0.6422751  0.7203086
#> 6     0 -4.07991970 0.3376486 -4.7416987 -3.4181407
#> 7     1 -4.20529807 0.6010068 -5.3832497 -3.0273464
#> 8     2 -4.06033981 0.7271590 -5.4855453 -2.6351344
#> 9     3 -4.07107236 0.6724557 -5.3890612 -2.7530835
#> 10    4 -3.58751419 0.6612129 -4.8834677 -2.2915606
#> 11    5 -3.92600168 0.1549628 -4.2297231 -3.6222803
# a by-law some towns repeal
d$law <- as.integer(on & !(i == 1 & d$year >= 2022))
stgdid(d$rate, d$town, d$year, trt = d$law, meth = "dcdh", reps = 49,
       seed = 1)$overall
#>      term        est        se         lo         hi n
#> 1  effect -3.6698465 0.4276516 -4.5080282 -2.8316647 9
#> 2 placebo -0.0742311 0.3569557 -0.7738514  0.6253892 9
```
