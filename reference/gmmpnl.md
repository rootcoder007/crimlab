# Difference and system GMM for dynamic panels

The generalized method of moments estimators of a dynamic panel with
unit effects: the difference estimator of Arellano and Bond (1991) and
the system estimator of Arellano and Bover (1995) and Blundell and Bond
(1998), with the instrument matrix built as Roodman (2009) describes it.
Andreatta et al. (2023, pp. 1372-1374, table 3) fit the logged rate of
waste crime in Italian provinces this way: one-step system GMM with
robust errors, time effects, the lagged rate and four enforcement and
organised-crime variables instrumented GMM-style by their second and
third lags, the AR(1) and AR(2) tests, the Hansen test, the Hansen test
of the differenced equations alone and the difference-in-Hansen test of
the level instruments, with no more instruments than provinces.

## Usage

``` r
gmmpnl(
  y,
  x = NULL,
  unit,
  time,
  ylag = 1,
  gmm = NULL,
  iv = NULL,
  z = NULL,
  lagmin = 2,
  lagmax = Inf,
  collapse = FALSE,
  type = "system",
  trans = "fd",
  steps = 1,
  robust = TRUE,
  wc = NULL,
  small = FALSE,
  hmat = "ab",
  cons = TRUE,
  tdum = TRUE,
  ar = 2,
  level = 0.95
)

# S3 method for class 'gmmpnl'
plot(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- y:

  The outcome of each row: finite numbers.

- x:

  `NULL` or the regressors: a numeric matrix or data frame, one row per
  row of `y`; the column names name them. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), a `gmmpnl`
  result.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

- ylag:

  The lags of `y` among the regressors: a whole number from 0 (default
  1, the first-order model of Arellano and Bond 1991, eq 1).

- gmm:

  `NULL` or the names of columns of `x` instrumented GMM-style, as the
  lagged outcome is.

- iv:

  `NULL` (every column of `x` not in `gmm`) or the names of the columns
  of `x` that instrument themselves, IV-style; `character(0)` for none.

- z:

  `NULL` or excluded IV-style instruments: a numeric matrix or data
  frame, one row per row of `y`.

- lagmin, lagmax:

  The lags of the GMM-style instruments: one value, or one per GMM-style
  variable (the lagged outcome first); whole numbers from 1, `lagmax`
  also `Inf`. The defaults 2 and `Inf` take every lag valid for an
  endogenous variable (Arellano and Bond 1991, sec. 2; Roodman 2009,
  sec. 3.2).

- collapse:

  `TRUE` for one instrument per lag rather than per lag and period
  (Roodman 2009, sec. 3.2); default `FALSE`, the full set of Arellano
  and Bond.

- type:

  `"system"` (default, Blundell and Bond 1998; the default of Roodman's
  xtabond2) or `"diff"` (Arellano and Bond 1991).

- trans:

  `"fd"` (first differences, default, Arellano and Bond 1991) or `"fod"`
  (forward orthogonal deviations, Arellano and Bover 1995).

- steps:

  1 (default; Andreatta et al. 2023 and the default of xtabond2) or 2
  (the two-step estimator).

- robust:

  `TRUE` (default) for the robust one-step variance and Windmeijer's
  correction of the two-step variance; `FALSE` for the conventional
  one-step variance and the uncorrected two-step variance.

- wc:

  `NULL` (default: as `robust`) or `TRUE` or `FALSE` for
  Windmeijer's (2005) correction of the two-step variance; it must agree
  with `robust`.

- small:

  `FALSE` (default) or `TRUE` for the small-sample variance, t
  statistics and an F test, as above.

- hmat:

  `"ab"` (default), `"dpd"` or `"id"`: `H` of the one-step weight, as
  above (Roodman 2009, secs. 3.3-3.4).

- cons:

  `TRUE` (default) for a constant in the level equation of system GMM;
  ignored by difference GMM, which removes it.

- tdum:

  `TRUE` (default) for time dummies, which Roodman (2009, sec. 3.3)
  advises so that the errors are not correlated across units.

- ar:

  The highest order of the Arellano-Bond tests: a whole number from 0
  (default 2, as xtabond2).

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

An object of class `gmmpnl`, a list: `coef` (`term`, `est`, `se`, `z`
(`t` with `small`), `p`, `lo`, `hi`); `tests` (`test`, `stat`, `df`,
`p`: AR(1) to AR(`ar`), Sargan, Hansen and, for system GMM, the Hansen
test without the level GMM-style instruments and the
difference-in-Hansen test, or their Sargan counterparts with the
conventional one-step variance); `wald`, one row (`test`, `"Wald"` or
`"F"`, `stat`, `df1`, `df2`, `p`); `fit`, one row: `n` (the rows of the
transformed equation for difference GMM, of the levels for system GMM,
as xtabond2 counts), `groups`, `inst`, `omin`, `oavg`, `omax` (rows per
group), `steps`, `k` and `sigma2`; `vcov`; `inst`, the instruments'
names; `dropped`, the time dummies dropped; `dummies`, those kept;
`note`; `type`; `trans`; `robust`; `small`; and `level`.

## Details

The model is \`y_it = a_1 y_i,t-1 + ... + a_L y_i,t-L + x_it b + d_t +
u_i

- v_it`with`L =
  ylag`. Periods are the distinct values of `time`in order: lags count positions in that order, and a unit without a row in a period has a gap there. A row enters the levels when the unit has the`ylag` periods before it; the transformed equation is the first difference of two such rows in consecutive periods (`trans
  =
  "fd"`), or the forward orthogonal deviation of a row from the mean of the unit's later rows, scaled by `sqrt(T
  / (T + 1))` and dated one period later (`trans = "fod"\`, Arellano and
  Bover 1995, sec. 3; Roodman 2009, eq 21), which loses no more than the
  last row of a unit with gaps.

GMM-style instruments (Roodman 2009, sec. 3.2) are made for the lagged
outcome (when `ylag > 0`, first) and for the columns of `x` named in
`gmm`: for the transformed equation dated `t`, the levels of the
variable at `t - l` for `l` from `lagmin` to `lagmax`, zero where the
unit has no row, one column per period and lag, or one per lag with
`collapse = TRUE`. For the level equation of system GMM, the difference
at `t - lagmin + 1` (one column per period, or one); the other lags are
redundant (Roodman 2009, sec. 3.4). Lags count from the variable: the
default `lagmin = 2` takes `y_i,t-2` and earlier for the equation in
`t`, the standard choice for an endogenous variable, and `lagmin = 1`
treats a column as predetermined. IV-style instruments (one column each,
transformed in the transformed equation and in levels in the level
equation, Roodman 2009, eq 24) are the columns of `x` in `iv`, the
columns of `z`, the time dummies and, in system GMM, the constant (which
is in the level equation only). Columns of instruments that are zero
throughout are dropped; time dummies collinear with the other regressors
in the stacked equations are dropped, and named in `dropped`.

The one-step weight is `(sum_i Z_i' H_i Z_i)^-1`. With `hmat = "ab"`
(the default), `H_i` is what the errors' variance would be with
independent, homoskedastic idiosyncratic errors and no unit effects in
the level equation: `M M'` for the transformed rows (2 on the diagonal
and -1 next to it for differences, the identity for deviations), `M`
between transformed and level rows and the identity for the levels
(Roodman 2009, eqs 22-23 and 26); `"dpd"` zeroes the blocks between the
equations, and `"id"` takes the identity throughout, as Blundell and
Bond (1998) did. The two-step weight is the inverse of
`sum_i Z_i' e_i e_i' Z_i` with the one-step residuals. A singular weight
matrix is replaced by its Moore-Penrose inverse, and the result says so.
With `robust = TRUE` (the default) the one-step variance is the robust
sandwich (Roodman 2009, eq 15) and the two-step variance is corrected as
Windmeijer (2005) shows (Roodman 2009, eq 18); with `robust = FALSE` the
one-step variance is the conventional `s2 (X'Z A Z'X)^-1`, eq 14 with
`Var(z e) = s2 Z'HZ`, `A = (Z'HZ)^-1` and `s2` as for Sargan's statistic
below, and the two-step variance is the uncorrected `(X'Z A_2 Z'X)^-1`
(Roodman 2009, p. 96). `wc` is the older name of the two-step choice:
left `NULL` it follows `robust`, and a value that contradicts `robust`,
or `wc` with one step, is refused. With `small = TRUE` the variance is
multiplied by `n / (n - k)` (`n` the rows reported in `fit`, `k` the
coefficients), and the coefficients have t statistics and the joint test
is an F test, both on the groups less one (less nothing without a
constant) degrees of freedom. Roodman (2009, p. 123) describes the
option of xtabond2 only as "small-sample corrections to the covariance
matrix estimate, resulting in t-test instead of z-test statistics for
the coefficients and an F test instead of a Wald chi2 test for overall
fit"; the degrees of freedom are those of his outputs (F(16, 140) for
140 groups without a constant, p. 113; F(12, 139) with one, p. 118), and
the factor `n / (n - k)` is this function's reading, as the paper gives
none.

The Arellano-Bond statistics (Arellano and Bond 1991, eqs 8-9) are those
of serial correlation of orders 1 to `ar` in the differenced residuals
of the estimate reported, also after deviations, with the variance's
terms taken from that estimate's residuals, weight and variance (Roodman
2009, sec. 3.5), and with the error variance taken as `s2 H` for the
conventional one-step variance (Roodman 2009, n. 19, the covariances of
independent homoskedastic errors through the transformation and the
differences, whatever `hmat`); AR(2) is the one that tests the
instruments. Hansen's statistic is that of the two-step estimate,
whichever step is reported, as robust overidentification tests need
(Roodman 2009, sec. 2.5); Sargan's is that of the one-step residuals,
`e'Z (Z'HZ)^-1 Z'e / s2` with `s2` the sum of the squared one-step
residuals of the transformed equation over the sum of their diagonal
elements of `H` (this function's estimate of the variance; Arellano and
Bond 1991, p. 282, write the statistic with an estimate of it). In
system GMM the Hansen test without the GMM-style instruments of the
level equation re-estimates the model with the rest, weighted by the
corresponding block of the full two-step moment variance (Roodman 2009,
n. 24), and the difference-in-Hansen statistic is the full statistic
less that one, on as many degrees of freedom as there are instruments
left out. With the conventional one-step variance Sargan's statistic is
the one reported for the instruments, and the difference test is the
difference-in-Sargan, the model re-estimated with the block of `Z'HZ`
(Roodman 2009, n. 24). `wald` is the joint test of the coefficients
other than the constant. Every test is on `m - k` degrees of freedom or
fewer, `m` the instruments and `k` the coefficients. When the
instruments outnumber the groups, `note` says so (Roodman 2009, sec.
2.6). Inference is normal unless `small = TRUE`. Andreatta et al.
describe their instruments both as "their second and third lags" (p.
1373) and as "Lag 2 and earlier lags" (table 3, note); `lagmin` and
`lagmax` take either reading.

The plot draws the coefficients other than the constant and the time
dummies, with their intervals, and returns what it drew.

## References

Arellano, M. and Bond, S. (1991). Some tests of specification for panel
data: Monte Carlo evidence and an application to employment equations.
*The Review of Economic Studies*, 58(2), 277-297.
[doi:10.2307/2297968](https://doi.org/10.2307/2297968)

Arellano, M. and Bover, O. (1995). Another look at the instrumental
variable estimation of error-components models. *Journal of
Econometrics*, 68(1), 29-51.
[doi:10.1016/0304-4076(94)01642-D](https://doi.org/10.1016/0304-4076%2894%2901642-D)

Blundell, R. and Bond, S. (1998). Initial conditions and moment
restrictions in dynamic panel data models. *Journal of Econometrics*,
87(1), 115-143.
[doi:10.1016/S0304-4076(98)00009-8](https://doi.org/10.1016/S0304-4076%2898%2900009-8)

Windmeijer, F. (2005). A finite sample correction for the variance of
linear efficient two-step GMM estimators. *Journal of Econometrics*,
126(1), 25-51.
[doi:10.1016/j.jeconom.2004.02.005](https://doi.org/10.1016/j.jeconom.2004.02.005)

Roodman, D. (2009). How to do xtabond2: an introduction to difference
and system GMM in Stata. *The Stata Journal*, 9(1), 86-136.
[doi:10.1177/1536867X0900900106](https://doi.org/10.1177/1536867X0900900106)

Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
Digging into waste: an analysis of waste crime in the Italian provinces.
*Regional Studies*, 57(7), 1367-1379.
[doi:10.1080/00343404.2022.2135697](https://doi.org/10.1080/00343404.2022.2135697)

## Examples

``` r
# 60 provinces over 8 years: a logged rate with persistence 0.5 and one
# endogenous driver, simulated
set.seed(7)
d <- expand.grid(year = 2010:2017, prov = paste0("p", 1:60))
u <- rnorm(60)[as.integer(substring(d$prov, 2))]
d$enf <- u + rnorm(480)
d$y <- 0
for (i in which(d$year > 2010)) {
  d$y[i] <- 0.5 * d$y[i - 1] + 0.3 * d$enf[i] + u[i] + rnorm(1)
}
r <- gmmpnl(d$y, d["enf"], d$prov, d$year, gmm = "enf", lagmax = 3,
            collapse = TRUE)
r
#> gmmpnl: system GMM (fd), step 1, 420 rows, 60 groups, 13 instruments
#> robust variance
#>     term         est         se           z            p         lo        hi
#>     L1.y  0.63565326 0.06952086  9.14334550 6.054957e-20  0.4993949 0.7719116
#>      enf  0.92497344 0.19159312  4.82780094 1.380490e-06  0.5494578 1.3004891
#>  (const) -0.06954309 0.13468780 -0.51632809 6.056253e-01 -0.3335263 0.1944401
#>    t2011  0.26161153 0.23231848  1.12609006 2.601274e-01 -0.1937243 0.7169474
#>    t2012 -0.01923081 0.19584762 -0.09819271 9.217793e-01 -0.4030851 0.3646235
#>    t2013 -0.03271092 0.19746408 -0.16565506 8.684284e-01 -0.4197334 0.3543116
#>    t2014  0.19530957 0.21204614  0.92107109 3.570133e-01 -0.2202932 0.6109124
#>    t2015  0.04273779 0.21878075  0.19534529 8.451226e-01 -0.3860646 0.4715402
#>    t2016  0.08698987 0.18530371  0.46944485 6.387517e-01 -0.2761987 0.4501785
#> Wald chi2(8) = 1275.2, p 5.394e-270
#>                      test       stat df            p
#>                     AR(1) -5.0629586 NA 4.127995e-07
#>                     AR(2)  0.6919913 NA 4.889428e-01
#>                    Sargan  8.8777291  4 6.422940e-02
#>                    Hansen  9.3189947  4 5.360236e-02
#>  Hansen without level GMM  0.6621633  2 7.181465e-01
#>      Difference-in-Hansen  8.6568314  2 1.318843e-02
r$tests
#>                       test       stat df            p
#> 1                    AR(1) -5.0629586 NA 4.127995e-07
#> 2                    AR(2)  0.6919913 NA 4.889428e-01
#> 3                   Sargan  8.8777291  4 6.422940e-02
#> 4                   Hansen  9.3189947  4 5.360236e-02
#> 5 Hansen without level GMM  0.6621633  2 7.181465e-01
#> 6     Difference-in-Hansen  8.6568314  2 1.318843e-02
```
