# Common correlated effects for panels with unobserved common factors

The common correlated effects estimators of Pesaran (2006) for a panel
whose errors hold unobserved common factors with unit-specific loadings,
correlated with the regressors: the dependent variable and the
regressors averaged across units in each period (with the weights `w`)
are added to each unit's regression, with the observed common factors
`d_t` (a constant, with `cons = TRUE`, and the columns of `fac`). With
`H = (D, ybar_w, Xbar_w)` and `M = I - H (H'H)^- H'` (eq 15, a
generalized inverse), the unit estimates are
`b_i = (X_i'M X_i)^-1 X_i'M y_i`. The mean group estimator `CCEMG` is
their mean (eq 53), with the nonparametric variance
`sum (b_i - b_MG)(b_i - b_MG)' / (N (N - 1))` (eq 58). The pooled
estimator `CCEP` is `(sum w_i X_i'M X_i)^-1 sum w_i X_i'M y_i`, the
pooling weights equal to the averaging weights as Theorem 3 requires,
with the variance of eqs 67-69, `sum w_i^2 Psi^-1 R Psi^-1`,
`Psi = sum w_i X_i'M X_i / T` and
`R = sum wt_i^2 (X_i'M X_i / T) (b_i - b_MG)(b_i - b_MG)' (X_i'M X_i / T) / (N - 1)`,
`wt_i = w_i / sqrt(mean(w^2))`. Inference is normal.

## Usage

``` r
ccepnl(
  y,
  x,
  unit,
  time,
  fac = NULL,
  cons = TRUE,
  w = NULL,
  avg = NULL,
  ylag = 0,
  lags = NULL,
  bias = "none",
  pool = NULL,
  pvar = "pesaran",
  level = 0.95
)

# S3 method for class 'ccepnl'
plot(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- y:

  The count of each row (a unit in a period): numbers from 0.

- x:

  The regressors: a numeric matrix or data frame, one row per row of
  `y`. For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), a
  `ccepnl` result.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

- fac:

  `NULL` or observed common factors: a numeric matrix or data frame, one
  row per row of `y`, constant within each period.

- cons:

  `TRUE` (default) for a constant among the observed factors,
  unit-specific intercepts.

- w:

  `NULL` (equal weights `1 / N`, as Pesaran 2006, sec. 5.2, says is
  typical) or a positive weight of each row, constant within each unit;
  they are scaled to sum to 1 over the units.

- avg:

  `NULL` (default: the averages of `y` and of every column of `x`,
  Pesaran's choice) or the averages kept among the columns of `H`: `"y"`
  and names of columns of `x`.

- ylag:

  The lags of `y` among the regressors: a whole number from 0 (default
  0, the static model of Pesaran 2006).

- lags:

  `NULL` (default: 0 for the static model, the integer part of `T^(1/3)`
  with `ylag > 0`) or the lags of the averages: a whole number from 0.

- bias:

  `"none"` (default), `"jack"` (half-panel jackknife) or `"rma"`
  (recursive mean adjustment): the correction of the mean group
  estimate, as above.

- pool:

  `NULL` (default: CCEP for the static model only), `TRUE` or `FALSE`:
  whether to make CCEP; `TRUE` is refused with `ylag > 0` or `bias`.

- pvar:

  `"pesaran"` (default, eqs 67-69 of Pesaran 2006) or `"bsw"` (Brown,
  Schmidt and Wooldridge 2026): the variance of CCEP.

- level:

  The confidence level of the intervals.

- ...:

  Passed to the base-graphics call.

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

An object of class `ccepnl`, a list: `mg` and `pooled`, the CCEMG and
CCEP estimates (`term`, `est`, `se`, `z`, `p`, `lo`, `hi`; `pooled`
`NULL` when not made); `mg0`, the uncorrected mean group estimate with
`bias` (`NULL` otherwise); `bi`, the unit estimates (units by
regressors, corrected with `bias`); `vmg`, `vp` and `vmg0`, the
variances; `hb`, the heterogeneity-bias test of CCEP (`NULL` without
CCEP); `nfac`, the number of factors (`er`, `thr`, `avg`); `fit`, one
row (`units`, `periods`, `k`, `h` the columns of `H`, `n` the rows of
the regressions, `lags`, `ylag`); `w`, the weights by unit; `bias`;
`pvar`; and `level`.

## Details

Pesaran (2006) derives the estimators for a balanced panel. An
unbalanced one is handled as standard implementations handle it (Ditzen
2018, p. 586, says only that xtdcce2 "supports unbalanced panels"; Su,
Wang and Wang 2026, p. 2, that "calculating cross-sectional averages of
the dependent variable and covariates remains straightforward in
unbalanced panels"): the averages of a period are taken over the units
observed in that period, their weights scaled to sum to 1 over those
units, and each unit's regression runs over its own periods, with its
rows of `H`. CCEMG and its variance are then as above. CCEP is
`(sum w_i X_i'M_i X_i)^-1 sum w_i X_i'M_i y_i`, and in its variance each
`X_i'M_i X_i` is divided by the unit's own number of periods, which is
this function's reading of `T` there. On a balanced panel all of this is
Pesaran's estimator. Su, Wang and Wang (2026, sec. 6.3, p. 20) compute
CCE on unbalanced panels in the same way, after Zhou and Zhang (2016, eq
12), and caution (p. 2 and n. 7, p. 4) that when whether a unit is
observed is related to its loadings and varies over time, the averages
no longer proxy the factors well; the estimators for unbalanced panels
with interactive effects of Czarnowske and Stammann (2026, Assumption 1
and Remark 1, p. 1309) likewise assume that the data are missing at
random, which excludes selection related to the idiosyncratic errors.
Neither paper gives a diagnostic or correction for CCE estimators, and
this function makes none. A period observed for one unit only, or a unit
with fewer periods than its regression needs (twice the regressors and
lags of `y`, the averages kept with their lags and the observed
factors), is refused by name.

**Dynamic panels.** With `ylag > 0` the lags `1..ylag` of `y` join the
regressors (as `L1.y`, ...), and with `lags = p` the averages kept enter
`H` with their lags `1..p` as well: the cross-sectionally augmented
regression of Chudik and Pesaran (2015, eq 24) for heterogeneous panels
with a lagged dependent variable and weakly exogenous regressors (lags
of `x`, as their eq 1 has, are columns of `x`). A unit's regression runs
over the periods in which it is observed together with its `ylag`
periods before, from the `p + 1`-th period on (eqs 25-26); the averages
are those of every observed unit-period, as above. `lags = NULL` takes
`p = 0` for the static model and, with `ylag > 0`, the integer part of
`T^(1/3)`, `T` the periods, as Chudik and Pesaran recommend (sec. 5.2,
eq 40, and p. 404: "we recommend setting `pT = [T^1/3]`"). The mean
group estimator is the mean of the unit estimates (eq 27) with the
nonparametric variance of eq 32 over `N`. Pooled dynamic CCE is not
supported: Chudik and Pesaran (2015, n. 5) recall that with
heterogeneous slopes "pooled estimators are inconsistent in the case of
panel data models with lagged dependent variables" (Pesaran and Smith
1995), and their theory is that of the mean group estimator;
`pool = NULL` (default) makes CCEP only for the static model, and
`pool = TRUE` with `ylag > 0` or with `bias` is refused.

`bias` corrects the mean group estimate for the small-T time-series bias
as Chudik and Pesaran (2015, sec. 4.1) define the corrections. `"jack"`
is the half-panel jackknife (sec. 4.1.1, the unnumbered equation after
eq 32), `2 b_MG - (b_MG^a + b_MG^b) / 2`, `a` and `b` the estimates over
periods `1..[T/2]` and `[T/2]+1..T`; each half is taken as a panel of
its own, so that the lags of `y` and of the averages of its first rows
come from inside the half, as "computed from the first half of the
available time period" reads, with the same `p`. `"rma"` is recursive
mean adjustment (sec. 4.1.2 and n. 13): `y` and each column of `x` less
the mean of the unit's values in the periods before,
`y_it - (t - 1)^-1 sum_{s<t} y_is`, a unit's first period dropped, and
the estimator applied to these series, so that the lags of `y` are lags
of the adjusted `y` and the averages are those of the adjusted values;
observed factors are not adjusted (Ditzen 2018, p. 589, leaves the
constant out). The paper states no variance for the corrected
estimators: `vmg` is the nonparametric mean group variance (eq 32)
applied to the corrected unit estimates (`bi`:
`2 b_i - (b_i^a + b_i^b) / 2` for the jackknife, the estimates on the
adjusted series for `"rma"`), and the uncorrected estimate and its
variance are returned in `mg0` and `vmg0` for comparison.

**Five lessons.** Juodis and Reese (2026) draw five lessons for applied
work. (1) The averages are factor estimates measured with error, which
can bias the estimates and change their variance: results should be
checked against the number and kind of averages (`avg` keeps a subset),
and with slopes common to all units and `T` not small against `N` the
intervals of CCEP are not centred without a bias correction, which this
function does not make (sec. 3.1). `nfac` reports the number of factors
that the averages kept identify, by the two regularized procedures of
their appendix A.1 on a balanced panel: `er`, the eigenvalue ratio of
Ahn and Horenstein (2013) on the averages normalised by `Sigma^-1/2`,
`Sigma = sum_i (Z_i - Zbar)'(Z_i - Zbar) / (N T)` (Juodis 2022; steps
1-3, `r_max` the averages less one); `thr`, the number of eigenvalues of
`Zbar'Zbar / T` at least 0.01 times the largest (Vogt et al. 2025; steps
1-2); and `avg`, the averages kept; `NA` where undefined. More averages
than factors is the case in which the paper warns of extra bias; the
regularized estimators themselves are not made. (2) A two-way
fixed-effects estimator
([`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md))
may still be consistent and should be compared and tested. `hb` is their
test of heterogeneity bias in CCEP (sec. 3.2, case study 2):
`delta_i = (X_i'M_i X_i - N^-1 sum_j X_j'M_j X_j) (b_i - b_MG)`,
`HB = delta' V^-1 delta` with `delta` their mean and
`V = N^-2 sum delta_i delta_i'`, chi-squared on `k` degrees of freedom
(`stat`, `df`, `p`); a rejection says that CCEP is biased by slopes that
vary with the regressors. (3) Dynamic models need the dynamics of every
observed variable specified, and a bias correction is recommended in
short panels (sec. 3.3; `ylag`, `lags`, `bias` above). (4) The CIPS
panel unit-root test does not test what CCE estimation needs; they
recommend the CSB test of Pesaran et al. (2013) (sec. 3.4), which this
package does not provide. (5) There is no general recipe for nonlinear
models, where the average of `y` should not be among the proxies (sec.
3.5); this function is linear. For cross-sectional dependence in the
residuals see
[`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md).

**Fixed-T variance.** Brown, Schmidt and Wooldridge (2026, eqs 12-14 and
Theorem 1) show that CCEP is the pooled least-squares estimate of `y` on
the regressors with the columns of `H` partialled out unit by unit, and
that time-varying regressors `d_t` added to it get coefficients of
exactly 0. Their fixed-T variance of CCEP, which counts the sampling
error of the averages (sec. 5.4, eqs 36-40, with `Psi = H` and
`q_i = vec(Z_i - Zbar)`, 0 in the columns of the observed factors), is
`pvar = "bsw"`: `A^-1 B A^-1 / N`, `A = N^-1 sum X_i'M X_i`,
`B = N^-1 sum s_i s_i'`,
`s_i = X_i'M u_i - G'(I + K)[Psi (Psi'Psi)^-1 (x) M] q_i`,
`u_i = y_i - X_i b`, `G = N^-1 sum u_i (x) (X_i - Xbar)` and `K` the
commutation matrix (eq 33); it is for the static pooled estimator with
equal weights on a balanced panel, and refused otherwise.

The plot draws the CCEMG and, when made, CCEP estimates of each
regressor with their intervals and returns what it drew.

## References

Pesaran, M. H. (2006). Estimation and inference in large heterogeneous
panels with a multifactor error structure. *Econometrica*, 74(4),
967-1012.
[doi:10.1111/j.1468-0262.2006.00692.x](https://doi.org/10.1111/j.1468-0262.2006.00692.x)

Ditzen, J. (2018). Estimating dynamic common-correlated effects in
Stata. *The Stata Journal*, 18(3), 585-617.
[doi:10.1177/1536867X1801800306](https://doi.org/10.1177/1536867X1801800306)

Su, L., Wang, F. and Wang, Y. (2026). Estimation and inference for
unbalanced panel data models with interactive fixed effects. *Journal of
Econometrics*, 255, 106222.
[doi:10.1016/j.jeconom.2026.106222](https://doi.org/10.1016/j.jeconom.2026.106222)

Juodis, A. and Reese, S. (2026). Five lessons for applied researchers
from twenty years of common correlated effects estimation. *Journal of
Econometrics*, 253, 106120.
[doi:10.1016/j.jeconom.2025.106120](https://doi.org/10.1016/j.jeconom.2025.106120)

Chudik, A. and Pesaran, M. H. (2015). Common correlated effects
estimation of heterogeneous dynamic panel data models with weakly
exogenous regressors. *Journal of Econometrics*, 188(2), 393-420.
[doi:10.1016/j.jeconom.2015.03.007](https://doi.org/10.1016/j.jeconom.2015.03.007)

Brown, N. L., Schmidt, P. and Wooldridge, J. M. (2026). A simple
reformulation of the common correlated effects model. *Econometric
Reviews*, 45(5), 760-772.
[doi:10.1080/07474938.2026.2615434](https://doi.org/10.1080/07474938.2026.2615434)

Czarnowske, D. and Stammann, A. (2026). Inference in unbalanced panel
data models with interactive fixed effects. *Econometric Reviews*,
45(9), 1306-1337.
[doi:10.1080/07474938.2026.2681519](https://doi.org/10.1080/07474938.2026.2681519)

## Examples

``` r
# 30 provinces over 12 years, one unobserved factor in y and x
set.seed(9)
d <- expand.grid(year = 1:12, prov = paste0("p", 1:30))
f <- rnorm(12)[d$year]
g <- rnorm(30, 1)[as.integer(substring(d$prov, 2))]
d$x <- g * f + rnorm(360)
d$y <- 1 + 0.7 * d$x + g * f + rnorm(360)
r <- ccepnl(d$y, d["x"], d$prov, d$year)
r
#> ccepnl: common correlated effects, 30 units, 12 periods
#> mean group (CCEMG):
#>  term       est         se        z            p        lo        hi
#>     x 0.7534011 0.06556476 11.49095 1.464985e-30 0.6248965 0.8819056
#> pooled (CCEP):
#>  term       est         se        z            p        lo        hi
#>     x 0.7581987 0.06456047 11.74401 7.580526e-32 0.6316625 0.8847349
#> heterogeneity bias of CCEP: chi2(1) = 0.0165, p 0.898
r$nfac
#>  er thr avg 
#>   1   2   2 
# a dynamic panel: y depends on its lag; jackknife-corrected mean group
set.seed(10)
tt <- 30
d <- expand.grid(year = 1:tt, prov = paste0("p", 1:40))
f <- rnorm(tt)[d$year]
g <- rnorm(40, 1)[as.integer(substring(d$prov, 2))]
d$x <- g * f + rnorm(nrow(d))
d$y <- 0
for (i in which(d$year > 1)) {
  d$y[i] <- 0.5 * d$y[i - 1] + 0.7 * d$x[i] + g[i] * f[i] + rnorm(1)
}
ccepnl(d$y, d["x"], d$prov, d$year, ylag = 1, bias = "jack")
#> Error: ccepnl: unit p1 has 12 usable periods in the first half of the periods: its regression on the regressors, the lags of y, the averages with their lags and the observed factors needs at least 13.
```
