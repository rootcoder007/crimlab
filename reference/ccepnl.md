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
  level = 0.95
)

# S3 method for class 'ccepnl'
plot(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- y:

  The outcome of each row: finite numbers.

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
CCEP estimates (`term`, `est`, `se`, `z`, `p`, `lo`, `hi`); `bi`, the
unit estimates (units by regressors); `vmg` and `vp`, the variances;
`fit`, one row (`units`, `periods`, `k`, `h` the columns of `H`, `n` the
rows); `w`, the weights by unit; and `level`.

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
Pesaran's estimator. Su, Wang and Wang (2026, p. 2) caution that when
whether a unit is observed is related to the regressors, the factors or
the loadings, the averages no longer proxy the factors well. A period
observed for one unit only, or a unit with fewer periods than its
regression needs (twice the regressors, the averages kept and the
observed factors), is refused by name.

Juodis and Reese (2026, sec. 1) draw two lessons this function leaves to
the user: the averages used can bias the estimates and their variance,
so results should be checked against the number and kind of averages
(`avg` keeps a subset; lesson 1), and under Pesaran's assumptions a
two-way fixed-effects estimator
([`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md))
may still be consistent and is worth comparing (lesson 2). With slopes
common to all units and `T` not small against `N`, the intervals of CCEP
are not centred without a bias correction, which this function does not
make (lesson 1, sec. 3.1). For cross-sectional dependence in the
residuals see
[`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md).

The plot draws the CCEMG and CCEP estimates of each regressor with their
intervals and returns what it drew.

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
```
