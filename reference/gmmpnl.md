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

  The count of each row (a unit in a period): numbers from 0.

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
  statistics and an F test as xtabond2 3.7.2 computes them (above).

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
`"F"`, `stat`, `df1`, `df2`, `p`); `xwald`, the same as xtabond2 prints
it (over every coefficient); `resid`, `y - x b` in levels, one per row
of `y`; `fit`, one row: `n` (the rows of the transformed equation for
difference GMM, of the levels for system GMM, as xtabond2 counts),
`groups`, `inst`, `omin`, `oavg`, `omax` (rows per group), `steps`, `k`
and `sigma2`; `vcov`; `inst`, the instruments' names; `dropped`, the
time dummies dropped; `dummies`, those kept; `note`; `type`; `trans`;
`robust`; `small`; and `level`.

## Details

The model is
`y_it = a_1 y_i,t-1 + ... + a_L y_i,t-L + x_it b + d_t + u_i + v_it`
with `L = ylag`. Periods are the distinct values of `time` in order:
lags count positions in that order, and a unit without a row in a period
has a gap there. A row enters the levels when the unit has the `ylag`
periods before it; the transformed equation is the first difference of
two such rows in consecutive periods (`trans = "fd"`), or the forward
orthogonal deviation of a row from the mean of the unit's later rows,
scaled by `sqrt(T / (T + 1))` and dated one period later
(`trans = "fod"`, Arellano and Bover 1995, sec. 3; Roodman 2009, eq 21),
which loses no more than the last row of a unit with gaps.

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
or `wc` with one step, is refused.

Roodman (2009, p. 123) describes `small` of xtabond2 only as
"small-sample corrections to the covariance matrix estimate, resulting
in t-test instead of z-test statistics for the coefficients and an F
test instead of a Wald chi2 test for overall fit"; `small = TRUE` makes
the corrections as xtabond2 3.7.2 computes them. With `n` the rows of
the equation reported (`fit$n`), `w` the rows of the transformed
equation in system GMM with `hmat` other than `"id"` and `n` otherwise,
`G` the groups, `k` the coefficients and `c` 1 for the constant of
system GMM (0 otherwise): the conventional one-step variance is
multiplied by `w / (w - k)`, and the robust one-step variance and both
two-step variances by `(n - 1) / (n - k) * G / (G - 1)`; `sigma2` is
multiplied by `w / (w - k)`. The coefficients have t statistics on
`df_r` degrees of freedom and the joint test is `F = W / q` on `q` and
`df_r`, with `df_r = n - q - c` for the conventional one-step variance
and `G - c` otherwise. The Sargan and Hansen statistics are unchanged;
so are the Arellano-Bond statistics, except after the conventional
one-step variance, where the program scales `sigma2` in the errors'
variance by `w / (w - k)` and leaves the coefficients' variance in the
statistic unscaled, and so does this function. Roodman's program also
lets the variance be clustered on other variables, which multiplies each
cluster's terms of the moments' variance by `G / (G - 1)` and drops that
factor from the variance; this function clusters by unit only.

`wald` is the joint test that the coefficients other than the constant
are 0, `W = b' V^-1 b` over them, on `q` degrees of freedom counted as
xtabond2 3.7.2 counts them, with or without `small`: the rank of the
regressors in the levels (system GMM; the transformed equation is left
out of the test) or in the transformed equation, less `c`; without a
constant, one less again when the vector of ones `i` lies in their
column space, which the program finds when `i'X_p (X'X)^-1 X_p'i`, `X_p`
those rows and `X` the whole stack, is within the machine epsilon times
the groups times the periods of the rows of `X_p`, relatively. The
program's own statistic takes `b' V^-1 b` over every coefficient, the
constant's included, while counting `q` without it; this function leaves
the constant out of the statistic too, so that the two agree when no
constant is estimated. `xwald` is the program's statistic, `b' V^-1 b`
over every coefficient on the same `q` (and `F = W / q` with `small`),
so that its printed line can be reproduced;
[`print()`](https://rdrr.io/r/base/print.html) shows it on a second line
when it differs from `wald`.

Two values differ from the program's by design. `sigma2` is that of the
one-step residuals, which Sargan's statistic uses, also with
`steps = 2`, where xtabond2 reports the two-step residuals' value. With
the conventional one-step variance the Arellano-Bond statistics take the
error variance as `s2` times the covariances of independent errors
through the differences whatever `hmat` (Roodman 2009, n. 19), where
with `hmat = "id"` in system GMM the program takes `sigma2` from the
levels and `H = I` there; the identity is the weight of Blundell and
Bond's first step, not the covariance of differenced errors.

`resid` holds `y - x b` in levels for every row of the level equation
(`NA` for rows without the `ylag` periods before), the unit effects
included; Andreatta et al. (2023, p. 1376) pass such residuals of the
dynamic model to Pesaran's test of cross-sectional dependence
([`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md))
and to Moran's I for each year
([`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)).

**Andreatta et al. (2023), step by step.** Their table 3 is one-step
system GMM with robust errors "with small-sample adjustments" (note):
`type = "system"`, `steps = 1`, `robust = TRUE`, `small = TRUE`; time
effects, `tdum = TRUE`; the lagged rate, `ylag = 1`, and the four
enforcement and organised-crime variables in `gmm`, instrumented by
"their second and third lags" (p. 1373), `lagmin = 2, lagmax = 3`, or by
"Lag 2 and earlier lags" (table note), `lagmax = Inf`, the other drivers
IV-style (`iv`); the difference of lag 1 in the level equation is what
`lagmin = 2` gives there. Their tests are the rows of `tests`: AR(1) and
AR(2), Hansen, "Hansen-level", the Hansen test of the first-differenced
model ("Hansen without level GMM"), and "Hansen-diff"
("Difference-in-Hansen"); "No. of instruments" is `fit$inst`, kept below
the 102 provinces (`note` warns otherwise). The supplementary "standard"
model with first lags as instruments is `lagmin = 1` for the variables
taken as predetermined. Before choosing a dynamic model they ran
Wooldridge's test on the static model's residuals
([`wsctst()`](https://rootcoder007.github.io/crimlab/reference/wsctst.md),
n. 9, F(1, 101)), and they compare the coefficient of the lagged rate
with pooled least squares ([`lm()`](https://rdrr.io/r/stats/lm.html))
and fixed effects
([`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md)
with the lag among the regressors). Their checks of spatial dependence
are
[`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)
and
[`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md)
on the rate with several weight matrices,
[`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md)
on the residuals of the static and dynamic models and
[`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)
on each year's residuals (p. 1376); their structural breaks (2013, 2014,
2015) are step dummies among the IV-style columns of `x`, and their
principal components of the enforcement rates are
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html).

Demir and Ozcan (2023, table 4) fit system GMM to 37 OECD countries over
2003-2017: the log rate of persons held in prison on its lag, the log
forest area per person and the log CO2 emissions per head, with the
controls the tests found endogenous instrumented GMM-style; where AR(2)
rejected, "the instruments were used at least with their third lags"
(`lagmin = 3`, with `ar = 3`), and the instruments were kept no more
than the countries (`collapse`, `lagmax`; `note` warns otherwise), with
the Hansen test above 0.10.

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
other than the constant (above). Every test is on `m - k` degrees of
freedom or fewer, `m` the instruments and `k` the coefficients. When the
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

Demir, C. and Ozcan, S. E. (2023). Do air quality and green space reduce
propensity to crime?: an empirical investigation of OECD countries.
*Management of Environmental Quality*, 34(2), 351-367.
[doi:10.1108/MEQ-08-2022-0235](https://doi.org/10.1108/MEQ-08-2022-0235)

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
#> Wald chi2(8) = 1560.9, p 0 (every coefficient, as xtabond2 prints it)
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

# the specification of Andreatta et al. (2023, table 3) on simulated
# data: the lagged rate and four drivers GMM-style by their second and
# third lags, collapsed so that the instruments stay below the
# provinces, two exogenous drivers, time effects, one-step robust
# errors with small-sample adjustments
set.seed(23)
np <- 102
d <- expand.grid(year = 2010:2019, prov = paste0("p", 1:np))
u <- rnorm(np)[as.integer(substring(d$prov, 2))]
n <- nrow(d)
for (v in c("clear", "conv", "mafia", "assoc")) d[[v]] <- u + rnorm(n)
d$gdp <- rnorm(n)
d$recyc <- rnorm(n)
d$wcrim <- 0
for (i in which(d$year > 2010)) {
  d$wcrim[i] <- 0.4 * d$wcrim[i - 1] - 0.1 * d$clear[i] +
    0.2 * d$mafia[i] + 0.1 * d$gdp[i] + u[i] + rnorm(1)
}
en <- c("clear", "conv", "mafia", "assoc")
a <- gmmpnl(d$wcrim, d[c(en, "gdp", "recyc")], d$prov, d$year, gmm = en,
            lagmin = 2, lagmax = 3, collapse = TRUE, small = TRUE)
a$tests
#>                       test      stat df            p
#> 1                    AR(1) -4.694829 NA 0.0000026683
#> 2                    AR(2)  1.432563 NA 0.1519828001
#> 3                   Sargan  9.750661 10 0.4626343573
#> 4                   Hansen 11.268943 10 0.3369524384
#> 5 Hansen without level GMM  1.868157  5 0.8670707506
#> 6     Difference-in-Hansen  9.400786  5 0.0941069902
a$fit$inst < a$fit$groups
#> [1] TRUE
# the static model's serial correlation and the dynamic residuals'
# cross-sectional dependence
wsctst(d$wcrim, d[c(en, "gdp", "recyc")], d$prov, d$year)
#> wsctst: Wooldridge test of serial correlation, F(1, 101) = 73.206, p 1.336e-13
#> lagged residual -0.2747 (se 0.02634; -0.5 without correlation), 816 pairs
k <- !is.na(a$resid)
pcdtst(a$resid[k], d$prov[k], d$year[k])
#>            cd             p           rho        absrho         pairs 
#> -2.063590e+00  3.905661e-02 -9.584211e-03  2.965807e-01  5.151000e+03 
#>         units       periods 
#>  1.020000e+02  9.000000e+00 
```
