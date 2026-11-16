#' Difference and system GMM for dynamic panels
#'
#' The generalized method of moments estimators of a dynamic panel with
#' unit effects: the difference estimator of Arellano and Bond (1991) and
#' the system estimator of Arellano and Bover (1995) and Blundell and Bond
#' (1998), with the instrument matrix built as Roodman (2009) describes
#' it. Andreatta et al. (2023, pp. 1372-1374, table 3) fit the logged rate
#' of waste crime in Italian provinces this way: one-step system GMM with
#' robust errors, time effects, the lagged rate and four enforcement and
#' organised-crime variables instrumented GMM-style by their second and
#' third lags, the AR(1) and AR(2) tests, the Hansen test, the Hansen test
#' of the differenced equations alone and the difference-in-Hansen test of
#' the level instruments, with no more instruments than provinces.
#'
#' The model is `y_it = a_1 y_i,t-1 + ... + a_L y_i,t-L + x_it b + d_t + u_i
#' + v_it` with `L = ylag`. Periods are the distinct values of `time` in
#' order: lags count positions in that order, and a unit without a row in
#' a period has a gap there. A row enters the levels when the unit has the
#' `ylag` periods before it; the transformed equation is the first
#' difference of two such rows in consecutive periods (`trans = "fd"`), or
#' the forward orthogonal deviation of a row from the mean of the unit's
#' later rows, scaled by `sqrt(T / (T + 1))` and dated one period later
#' (`trans = "fod"`, Arellano and Bover 1995, sec. 3; Roodman 2009, eq 21),
#' which loses no more than the last row of a unit with gaps.
#'
#' GMM-style instruments (Roodman 2009, sec. 3.2) are made for the lagged
#' outcome (when `ylag > 0`, first) and for the columns of `x` named in
#' `gmm`: for the transformed equation dated `t`, the levels of the
#' variable at `t - l` for `l` from `lagmin` to `lagmax`, zero where the
#' unit has no row, one column per period and lag, or one per lag with
#' `collapse = TRUE`. For the level equation of system GMM, the difference
#' at `t - lagmin + 1` (one column per period, or one); the other lags are
#' redundant (Roodman 2009, sec. 3.4). Lags count from the variable: the
#' default `lagmin = 2` takes `y_i,t-2` and earlier for the equation in
#' `t`, the standard choice for an endogenous variable, and `lagmin = 1`
#' treats a column as predetermined. IV-style instruments (one column each,
#' transformed in the transformed equation and in levels in the level
#' equation, Roodman 2009, eq 24) are the columns of `x` in `iv`, the
#' columns of `z`, the time dummies and, in system GMM, the constant
#' (which is in the level equation only). Columns of instruments that are
#' zero throughout are dropped; time dummies collinear with the other
#' regressors in the stacked equations are dropped, and named in
#' `dropped`.
#'
#' The one-step weight is `(sum_i Z_i' H_i Z_i)^-1`. With `hmat = "ab"`
#' (the default), `H_i` is what the errors' variance would be with
#' independent, homoskedastic idiosyncratic errors and no unit effects in
#' the level equation: `M M'` for the transformed rows (2 on the diagonal
#' and -1 next to it for differences, the identity for deviations), `M`
#' between transformed and level rows and the identity for the levels
#' (Roodman 2009, eqs 22-23 and 26); `"dpd"` zeroes the blocks between the
#' equations, and `"id"` takes the identity throughout, as Blundell and
#' Bond (1998) did. The two-step weight is the inverse of `sum_i Z_i' e_i
#' e_i' Z_i` with the one-step residuals. A singular weight matrix is
#' replaced by its Moore-Penrose inverse, and the result says so. With
#' `robust = TRUE` (the default) the one-step variance is the robust
#' sandwich (Roodman 2009, eq 15) and the two-step variance is corrected as
#' Windmeijer (2005) shows (Roodman 2009, eq 18); with `robust = FALSE` the
#' one-step variance is the conventional `s2 (X'Z A Z'X)^-1`, eq 14 with
#' `Var(z e) = s2 Z'HZ`, `A = (Z'HZ)^-1` and `s2` as for Sargan's statistic
#' below, and the two-step variance is the uncorrected `(X'Z A_2 Z'X)^-1`
#' (Roodman 2009, p. 96). `wc` is the older name of the two-step choice:
#' left `NULL` it follows `robust`, and a value that contradicts `robust`,
#' or `wc` with one step, is refused.
#'
#' Roodman (2009, p. 123) describes `small` of xtabond2 only as
#' "small-sample corrections to the covariance matrix estimate, resulting
#' in t-test instead of z-test statistics for the coefficients and an F
#' test instead of a Wald chi2 test for overall fit"; `small = TRUE` makes
#' the corrections as xtabond2 3.7.2 computes them. With `n` the rows of
#' the equation reported (`fit$n`), `w` the rows of the transformed
#' equation in system GMM with `hmat` other than `"id"` and `n` otherwise,
#' `G` the groups, `k` the coefficients and `c` 1 for the constant of
#' system GMM (0 otherwise): the conventional one-step variance is
#' multiplied by `w / (w - k)`, and the robust one-step variance and both
#' two-step variances by `(n - 1) / (n - k) * G / (G - 1)`; `sigma2` is
#' multiplied by `w / (w - k)`. The coefficients have t statistics on
#' `df_r` degrees of freedom and the joint test is `F = W / q` on `q` and
#' `df_r`, with `df_r = n - q - c` for the conventional one-step variance
#' and `G - c` otherwise. The Sargan and Hansen statistics are unchanged;
#' so are the Arellano-Bond statistics, except after the conventional one-step
#' variance, where the program scales `sigma2` in the errors' variance by
#' `w / (w - k)` and leaves the coefficients' variance in the statistic
#' unscaled, and so does this function. Roodman's program also lets the
#' variance be clustered on other variables, which multiplies each
#' cluster's terms of the moments' variance by `G / (G - 1)` and drops
#' that factor from the variance; this function clusters by unit only.
#'
#' `wald` is the joint test that the coefficients other than the constant
#' are 0, `W = b' V^-1 b` over them, on `q` degrees of freedom counted as
#' xtabond2 3.7.2 counts them, with or without `small`: the rank of the
#' regressors in the levels (system GMM; the transformed equation is left
#' out of the test) or in the transformed equation, less `c`; without a
#' constant, one less again when the vector of ones `i` lies in their
#' column space, which the program finds when `i'X_p (X'X)^-1 X_p'i`, `X_p`
#' those rows and `X` the whole stack, is within the machine epsilon times
#' the groups times the periods of the rows of `X_p`, relatively. The
#' program's own statistic takes `b' V^-1 b` over every coefficient,
#' the constant's included, while counting `q` without it; this function
#' leaves the constant out of the statistic too, so that the two agree
#' when no constant is estimated. `xwald` is the program's statistic, `b'
#' V^-1 b` over every coefficient on the same `q` (and `F = W / q` with
#' `small`), so that its printed line can be reproduced; `print()` shows
#' it on a second line when it differs from `wald`.
#'
#' Two values differ from the program's by design. `sigma2` is that of
#' the one-step residuals, which Sargan's statistic uses, also with
#' `steps = 2`, where xtabond2 reports the two-step residuals' value. With
#' the conventional one-step variance the Arellano-Bond statistics take
#' the error variance as `s2` times the covariances of independent errors
#' through the differences whatever `hmat` (Roodman 2009, n. 19), where
#' with `hmat = "id"` in system GMM the program takes `sigma2` from the
#' levels and `H = I` there; the identity is the weight of Blundell and
#' Bond's first step, not the covariance of differenced errors.
#'
#' `resid` holds `y - x b` in levels for every row of the level equation
#' (`NA` for rows without the `ylag` periods before), the unit effects
#' included; Andreatta et al. (2023, p. 1376) pass such residuals of the
#' dynamic model to Pesaran's test of cross-sectional dependence
#' ([pcdtst()]) and to Moran's I for each year ([gmoran()]).
#'
#' **Andreatta et al. (2023), step by step.** Their table 3 is one-step
#' system GMM with robust errors "with small-sample adjustments" (note):
#' `type = "system"`, `steps = 1`, `robust = TRUE`, `small = TRUE`; time
#' effects, `tdum = TRUE`; the lagged rate, `ylag = 1`, and the four
#' enforcement and organised-crime variables in `gmm`, instrumented by
#' "their second and third lags" (p. 1373), `lagmin = 2, lagmax = 3`, or
#' by "Lag 2 and earlier lags" (table note), `lagmax = Inf`, the other
#' drivers IV-style (`iv`); the difference of lag 1 in the level equation
#' is what `lagmin = 2` gives there. Their tests are the rows of `tests`:
#' AR(1) and AR(2), Hansen, "Hansen-level", the Hansen test of the
#' first-differenced model ("Hansen without level GMM"), and "Hansen-diff"
#' ("Difference-in-Hansen"); "No. of instruments" is `fit$inst`, kept
#' below the 102 provinces (`note` warns otherwise). The supplementary
#' "standard" model with first lags as instruments is `lagmin = 1` for
#' the variables taken as predetermined. Before choosing a dynamic model
#' they ran Wooldridge's test on the static model's residuals ([wsctst()],
#' n. 9, F(1, 101)), and they compare the coefficient of the lagged rate
#' with pooled least squares (`lm()`) and fixed effects ([fepnl()] with the
#' lag among the regressors). Their checks of spatial dependence are
#' [gmoran()] and [gearc()] on the rate with several weight matrices,
#' [pcdtst()] on the residuals of the static and dynamic models and
#' [gmoran()] on each year's residuals (p. 1376); their structural breaks
#' (2013, 2014, 2015) are step dummies among the IV-style columns of `x`,
#' and their principal components of the enforcement rates are
#' `stats::prcomp()`.
#'
#' The Arellano-Bond statistics (Arellano and Bond 1991, eqs 8-9) are
#' those of serial correlation of orders 1 to `ar` in the differenced
#' residuals of the estimate reported, also after deviations, with the
#' variance's terms taken from that estimate's residuals, weight and
#' variance (Roodman 2009, sec. 3.5), and with the error variance taken as
#' `s2 H` for the conventional one-step variance (Roodman 2009, n. 19, the
#' covariances of independent homoskedastic errors through the
#' transformation and the differences, whatever `hmat`); AR(2) is the one
#' that tests the instruments. Hansen's statistic is that of the two-step
#' estimate, whichever step is reported, as robust overidentification
#' tests need (Roodman 2009, sec. 2.5); Sargan's is that of the one-step
#' residuals, `e'Z (Z'HZ)^-1 Z'e / s2` with `s2` the sum of the squared one-step
#' residuals of the transformed equation over the sum of their diagonal
#' elements of `H` (this function's estimate of the variance; Arellano and
#' Bond 1991, p. 282, write the statistic with an estimate of it). In
#' system GMM the Hansen test without the GMM-style instruments of the level
#' equation re-estimates the model with the rest, weighted by the
#' corresponding block of the full two-step moment variance (Roodman 2009,
#' n. 24), and the difference-in-Hansen statistic is the full statistic
#' less that one, on as many degrees of freedom as there are instruments
#' left out. With the conventional one-step variance Sargan's statistic is
#' the one reported for the instruments, and the difference test is the
#' difference-in-Sargan, the model re-estimated with the block of `Z'HZ`
#' (Roodman 2009, n. 24). `wald` is the joint test of the coefficients
#' other than the constant (above). Every test is on `m - k` degrees of
#' freedom or fewer, `m` the instruments and `k` the coefficients. When the
#' instruments outnumber the groups, `note` says so (Roodman 2009, sec.
#' 2.6). Inference is normal unless `small = TRUE`. Andreatta et al.
#' describe their instruments both as "their second and third lags" (p.
#' 1373) and as "Lag 2 and earlier lags" (table 3, note); `lagmin` and
#' `lagmax` take either reading.
#'
#' @inheritParams lmdix
#' @inheritParams fepnl
#' @param x `NULL` or the regressors: a numeric matrix or data frame, one
#'   row per row of `y`; the column names name them. For `plot()`, a
#'   `gmmpnl` result.
#' @param ylag The lags of `y` among the regressors: a whole number from 0
#'   (default 1, the first-order model of Arellano and Bond 1991, eq 1).
#' @param gmm `NULL` or the names of columns of `x` instrumented
#'   GMM-style, as the lagged outcome is.
#' @param iv `NULL` (every column of `x` not in `gmm`) or the names of the
#'   columns of `x` that instrument themselves, IV-style; `character(0)`
#'   for none.
#' @param z `NULL` or excluded IV-style instruments: a numeric matrix or
#'   data frame, one row per row of `y`.
#' @param lagmin,lagmax The lags of the GMM-style instruments: one value,
#'   or one per GMM-style variable (the lagged outcome first); whole
#'   numbers from 1, `lagmax` also `Inf`. The defaults 2 and `Inf` take
#'   every lag valid for an endogenous variable (Arellano and Bond 1991,
#'   sec. 2; Roodman 2009, sec. 3.2).
#' @param collapse `TRUE` for one instrument per lag rather than per lag and
#'   period (Roodman 2009, sec. 3.2); default `FALSE`, the full set of
#'   Arellano and Bond.
#' @param type `"system"` (default, Blundell and Bond 1998; the default of
#'   Roodman's xtabond2) or `"diff"` (Arellano and Bond 1991).
#' @param trans `"fd"` (first differences, default, Arellano and Bond
#'   1991) or `"fod"` (forward orthogonal deviations, Arellano and Bover
#'   1995).
#' @param steps 1 (default; Andreatta et al. 2023 and the default of
#'   xtabond2) or 2 (the two-step estimator).
#' @param robust `TRUE` (default) for the robust one-step variance and
#'   Windmeijer's correction of the two-step variance; `FALSE` for the
#'   conventional one-step variance and the uncorrected two-step variance.
#' @param wc `NULL` (default: as `robust`) or `TRUE` or `FALSE` for
#'   Windmeijer's (2005) correction of the two-step variance; it must agree
#'   with `robust`.
#' @param small `FALSE` (default) or `TRUE` for the small-sample variance,
#'   t statistics and an F test as xtabond2 3.7.2 computes them (above).
#' @param hmat `"ab"` (default), `"dpd"` or `"id"`: `H` of the one-step
#'   weight, as above (Roodman 2009, secs. 3.3-3.4).
#' @param cons `TRUE` (default) for a constant in the level equation of
#'   system GMM; ignored by difference GMM, which removes it.
#' @param tdum `TRUE` (default) for time dummies, which Roodman (2009, sec.
#'   3.3) advises so that the errors are not correlated across units.
#' @param ar The highest order of the Arellano-Bond tests: a whole number
#'   from 0 (default 2, as xtabond2).
#' @return An object of class `gmmpnl`, a list: `coef` (`term`, `est`, `se`,
#'   `z` (`t` with `small`), `p`, `lo`, `hi`); `tests` (`test`, `stat`,
#'   `df`, `p`: AR(1) to AR(`ar`), Sargan, Hansen and, for system GMM, the
#'   Hansen test without the level GMM-style instruments and the
#'   difference-in-Hansen test, or their Sargan counterparts with the
#'   conventional one-step variance); `wald`, one row (`test`, `"Wald"` or
#'   `"F"`, `stat`, `df1`, `df2`, `p`); `xwald`, the same as xtabond2
#'   prints it (over every coefficient); `resid`, `y - x b` in levels, one
#'   per row of `y`;
#'   `fit`, one row: `n` (the rows of the transformed equation for
#'   difference GMM, of the levels for system GMM, as xtabond2 counts),
#'   `groups`, `inst`, `omin`, `oavg`, `omax` (rows per group), `steps`,
#'   `k` and `sigma2`; `vcov`; `inst`, the instruments' names; `dropped`,
#'   the time dummies dropped; `dummies`, those kept; `note`; `type`;
#'   `trans`; `robust`; `small`; and `level`.
#' @references
#'   Arellano, M. and Bond, S. (1991). Some tests of specification for
#'   panel data: Monte Carlo evidence and an application to employment
#'   equations. *The Review of Economic Studies*, 58(2), 277-297.
#'   \doi{10.2307/2297968}
#'
#'   Arellano, M. and Bover, O. (1995). Another look at the instrumental
#'   variable estimation of error-components models. *Journal of
#'   Econometrics*, 68(1), 29-51. \doi{10.1016/0304-4076(94)01642-D}
#'
#'   Blundell, R. and Bond, S. (1998). Initial conditions and moment
#'   restrictions in dynamic panel data models. *Journal of Econometrics*,
#'   87(1), 115-143. \doi{10.1016/S0304-4076(98)00009-8}
#'
#'   Windmeijer, F. (2005). A finite sample correction for the variance of
#'   linear efficient two-step GMM estimators. *Journal of Econometrics*,
#'   126(1), 25-51. \doi{10.1016/j.jeconom.2004.02.005}
#'
#'   Roodman, D. (2009). How to do xtabond2: an introduction to difference
#'   and system GMM in Stata. *The Stata Journal*, 9(1), 86-136.
#'   \doi{10.1177/1536867X0900900106}
#'
#'   Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
#'   Digging into waste: an analysis of waste crime in the Italian
#'   provinces. *Regional Studies*, 57(7), 1367-1379.
#'   \doi{10.1080/00343404.2022.2135697}
#' @examples
#' # 60 provinces over 8 years: a logged rate with persistence 0.5 and one
#' # endogenous driver, simulated
#' set.seed(7)
#' d <- expand.grid(year = 2010:2017, prov = paste0("p", 1:60))
#' u <- rnorm(60)[as.integer(substring(d$prov, 2))]
#' d$enf <- u + rnorm(480)
#' d$y <- 0
#' for (i in which(d$year > 2010)) {
#'   d$y[i] <- 0.5 * d$y[i - 1] + 0.3 * d$enf[i] + u[i] + rnorm(1)
#' }
#' r <- gmmpnl(d$y, d["enf"], d$prov, d$year, gmm = "enf", lagmax = 3,
#'             collapse = TRUE)
#' r
#' r$tests
#'
#' # the specification of Andreatta et al. (2023, table 3) on simulated
#' # data: the lagged rate and four drivers GMM-style by their second and
#' # third lags, collapsed so that the instruments stay below the
#' # provinces, two exogenous drivers, time effects, one-step robust
#' # errors with small-sample adjustments
#' set.seed(23)
#' np <- 102
#' d <- expand.grid(year = 2010:2019, prov = paste0("p", 1:np))
#' u <- rnorm(np)[as.integer(substring(d$prov, 2))]
#' n <- nrow(d)
#' for (v in c("clear", "conv", "mafia", "assoc")) d[[v]] <- u + rnorm(n)
#' d$gdp <- rnorm(n)
#' d$recyc <- rnorm(n)
#' d$wcrim <- 0
#' for (i in which(d$year > 2010)) {
#'   d$wcrim[i] <- 0.4 * d$wcrim[i - 1] - 0.1 * d$clear[i] +
#'     0.2 * d$mafia[i] + 0.1 * d$gdp[i] + u[i] + rnorm(1)
#' }
#' en <- c("clear", "conv", "mafia", "assoc")
#' a <- gmmpnl(d$wcrim, d[c(en, "gdp", "recyc")], d$prov, d$year, gmm = en,
#'             lagmin = 2, lagmax = 3, collapse = TRUE, small = TRUE)
#' a$tests
#' a$fit$inst < a$fit$groups
#' # the static model's serial correlation and the dynamic residuals'
#' # cross-sectional dependence
#' wsctst(d$wcrim, d[c(en, "gdp", "recyc")], d$prov, d$year)
#' k <- !is.na(a$resid)
#' pcdtst(a$resid[k], d$prov[k], d$year[k])
#' @export
gmmpnl <- function(y, x = NULL, unit, time, ylag = 1, gmm = NULL, iv = NULL,
                   z = NULL, lagmin = 2, lagmax = Inf, collapse = FALSE,
                   type = "system", trans = "fd", steps = 1, robust = TRUE,
                   wc = NULL, small = FALSE, hmat = "ab", cons = TRUE,
                   tdum = TRUE, ar = 2, level = 0.95) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(z)) z <- as.matrix(z)
  .Call(C_clab_r_gmmpnl, y, x, unit, time, ylag, gmm, iv, z, lagmin, lagmax,
        collapse, type, trans, steps, robust, wc, small, hmat, cons, tdum, ar,
        level)
}

#' @export
print.gmmpnl <- function(x, ...) {
  f <- x[["fit"]]
  cat("gmmpnl: ", if (x[["type"]] == "system") "system" else "difference",
      " GMM (", x[["trans"]], "), step ", f[["steps"]], ", ", f[["n"]],
      " rows, ", f[["groups"]], " groups, ", f[["inst"]], " instruments\n",
      sep = "")
  cat(if (isTRUE(x[["robust"]])) "robust" else "conventional",
      " variance", if (isTRUE(x[["small"]])) ", small-sample" else "", "\n",
      sep = "")
  print(x[["coef"]], row.names = FALSE)
  w <- x[["wald"]]
  if (!is.null(w)) cat(.gwline(w, ""))
  xw <- x[["xwald"]]
  if (!is.null(xw) && !identical(xw$stat, w$stat)) {
    cat(.gwline(xw, " (every coefficient, as xtabond2 prints it)"))
  }
  print(x[["tests"]], row.names = FALSE)
  for (s in x[["note"]]) cat("note: ", s, "\n", sep = "")
  invisible(x)
}

#' @rdname gmmpnl
#' @inheritParams plotlrnz
#' @details The plot draws the coefficients other than the constant and
#'   the time dummies, with their intervals, and returns what it drew.
#' @export
plot.gmmpnl <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "gmmpnl", "coef")
  r <- x$coef
  r <- r[!r$term %in% c("(const)", x$dummies), , drop = FALSE]
  .pdots(r$term, r$est, r$lo, r$hi, .pmain(main, "Dynamic panel GMM"),
         "Coefficient", col, pal, ref = 0, ...)
}

#' Wooldridge's test for serial correlation in a panel
#'
#' Whether the idiosyncratic errors of a linear panel model with unit
#' effects are serially correlated, by the test Wooldridge (2002) derives
#' and Drukker (2003) sets out and evaluates, the test Andreatta et al.
#' (2023, p. 1373) run on the residuals of their static models before
#' choosing a dynamic one. The first differences of `y` are regressed on
#' those of `x`, without a constant, which differencing removes (Drukker
#' 2003, sec. 2); without serial correlation the residuals `e_it` of
#' consecutive differences correlate at -0.5, so `e_it` is regressed on
#' `e_i,t-1` with the variance clustered by unit, and the coefficient is
#' tested against -0.5: `F = ((coef + 0.5) / se)^2` on 1 and `G - 1`
#' degrees of freedom, `G` the units with a pair. A difference needs the
#' period before, so a gap removes two differences and the pairs that use
#' them. Drukker does not say whether the second regression has a
#' constant; `cons = FALSE` (the default) leaves it out, since the
#' residuals of the first have none to absorb.
#'
#' @inheritParams gmmpnl
#' @inheritParams fepnl
#' @param x `NULL` or the regressors: a numeric matrix or data frame, one
#'   row per row of `y`.
#' @param cons `TRUE` for a constant in the regression of the residual on
#'   its lag.
#' @return An object of class `wsctst`, a list: `stat`, a named vector
#'   (`F`, `df1`, `df2`, `p`, `coef` the coefficient of the lagged residual,
#'   `se`, `pairs`, `clusters`, `ndiff` the differences); `coef`, the
#'   differenced regression (`term`, `est`, `se`), `NULL` without `x`; and
#'   `ssc`.
#' @references
#'   Drukker, D. M. (2003). Testing for serial correlation in linear
#'   panel-data models. *The Stata Journal*, 3(2), 168-177.
#'   \doi{10.1177/1536867X0300300206}
#'
#'   Wooldridge, J. M. (2002). *Econometric Analysis of Cross Section and
#'   Panel Data*. MIT Press.
#'
#'   Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
#'   Digging into waste: an analysis of waste crime in the Italian
#'   provinces. *Regional Studies*, 57(7), 1367-1379.
#'   \doi{10.1080/00343404.2022.2135697}
#' @examples
#' set.seed(8)
#' d <- expand.grid(year = 1:8, prov = paste0("p", 1:40))
#' d$x <- rnorm(320)
#' d$e <- 0
#' for (i in which(d$year > 1)) d$e[i] <- 0.6 * d$e[i - 1] + rnorm(1)
#' d$y <- d$x + rnorm(40)[as.integer(substring(d$prov, 2))] + d$e
#' wsctst(d$y, d["x"], d$prov, d$year)
#' @export
wsctst <- function(y, x = NULL, unit, time, cons = FALSE, ssc = "cr1") {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_wsctst, y, x, unit, time, cons, ssc)
}

#' @export
print.wsctst <- function(x, ...) {
  s <- x[["stat"]]
  cat("wsctst: Wooldridge test of serial correlation, F(1, ", s[["df2"]],
      ") = ", format(s[["F"]], digits = 5), ", p ",
      format(s[["p"]], digits = 4), "\n", sep = "")
  cat("lagged residual ", format(s[["coef"]], digits = 4), " (se ",
      format(s[["se"]], digits = 4), "; -0.5 without correlation), ",
      s[["pairs"]], " pairs\n", sep = "")
  invisible(x)
}

#' Common correlated effects for panels with unobserved common factors
#'
#' The common correlated effects estimators of Pesaran (2006) for a panel
#' whose errors hold unobserved common factors with unit-specific
#' loadings, correlated with the regressors: the dependent variable and
#' the regressors averaged across units in each period (with the weights
#' `w`) are added to each unit's regression, with the observed common
#' factors `d_t` (a constant, with `cons = TRUE`, and the columns of `fac`).
#' With `H = (D, ybar_w, Xbar_w)` and `M = I - H (H'H)^- H'` (eq 15, a
#' generalized inverse), the unit estimates are `b_i = (X_i'M X_i)^-1
#' X_i'M y_i`. The mean group estimator `CCEMG` is their mean (eq 53),
#' with the nonparametric variance `sum (b_i - b_MG)(b_i - b_MG)' / (N (N -
#' 1))` (eq 58). The pooled estimator `CCEP` is `(sum w_i X_i'M X_i)^-1 sum
#' w_i X_i'M y_i`, the pooling weights equal to the averaging weights as
#' Theorem 3 requires, with the variance of eqs 67-69, `sum w_i^2 Psi^-1 R
#' Psi^-1`, `Psi = sum w_i X_i'M X_i / T` and `R = sum wt_i^2 (X_i'M X_i /
#' T) (b_i - b_MG)(b_i - b_MG)' (X_i'M X_i / T) / (N - 1)`, `wt_i = w_i /
#' sqrt(mean(w^2))`. Inference is normal.
#'
#' Pesaran (2006) derives the estimators for a balanced panel. An
#' unbalanced one is handled as standard implementations handle it (Ditzen
#' 2018, p. 586, says only that xtdcce2 "supports unbalanced panels";
#' Su, Wang and Wang 2026, p. 2, that "calculating cross-sectional
#' averages of the dependent variable and covariates remains
#' straightforward in unbalanced panels"): the averages of a period are
#' taken over the units observed in that period, their weights scaled to
#' sum to 1 over those units, and each unit's regression runs over its own
#' periods, with its rows of `H`. CCEMG and its variance are then as above.
#' CCEP is `(sum w_i X_i'M_i X_i)^-1 sum w_i X_i'M_i y_i`, and in its
#' variance each `X_i'M_i X_i` is divided by the unit's own number of
#' periods, which is this function's reading of `T` there. On a balanced
#' panel all of this is Pesaran's estimator.
#' Su, Wang and Wang (2026, sec. 6.3, p. 20) compute CCE on unbalanced
#' panels in the same way, after Zhou and Zhang (2016, eq 12), and caution
#' (p. 2 and n. 7, p. 4) that when whether a unit is observed is related
#' to its loadings and varies over time, the averages no longer proxy the
#' factors well; the estimators for unbalanced panels with interactive
#' effects of Czarnowske and Stammann (2026, Assumption 1 and Remark 1, p.
#' 1309) likewise assume that the data are missing at random, which
#' excludes selection related to the idiosyncratic errors. Neither paper
#' gives a diagnostic or correction for CCE estimators, and this function
#' makes none. A period observed for one unit only, or a unit with fewer
#' periods than its regression needs (twice the regressors and lags of
#' `y`, the averages kept with their lags and the observed factors), is
#' refused by name.
#'
#' **Dynamic panels.** With `ylag > 0` the lags `1..ylag` of `y` join the
#' regressors (as `L1.y`, ...), and with `lags = p` the averages kept enter
#' `H` with their lags `1..p` as well: the cross-sectionally augmented
#' regression of Chudik and Pesaran (2015, eq 24) for heterogeneous panels
#' with a lagged dependent variable and weakly exogenous regressors (lags
#' of `x`, as their eq 1 has, are columns of `x`). A unit's regression runs
#' over the periods in which it is observed together with its `ylag`
#' periods before, from the `p + 1`-th period on (eqs 25-26); the averages
#' are those of every observed unit-period, as above. `lags = NULL` takes
#' `p = 0` for the static model and, with `ylag > 0`, the integer part of
#' `T^(1/3)`, `T` the periods, as Chudik and Pesaran recommend (sec. 5.2,
#' eq 40, and p. 404: "we recommend setting `pT = [T^1/3]`"). The mean group
#' estimator is the mean of the unit estimates (eq 27) with the
#' nonparametric variance of eq 32 over `N`. Pooled dynamic CCE is not
#' supported: Chudik and Pesaran (2015, n. 5) recall that with
#' heterogeneous slopes "pooled estimators are inconsistent in the case
#' of panel data models with lagged dependent variables" (Pesaran and
#' Smith 1995), and their
#' theory is that of the mean group estimator; `pool = NULL` (default)
#' makes CCEP only for the static model, and `pool = TRUE` with `ylag > 0`
#' or with `bias` is refused.
#'
#' `bias` corrects the mean group estimate for the small-T time-series
#' bias as Chudik and Pesaran (2015, sec. 4.1) define the corrections.
#' `"jack"` is the half-panel jackknife (sec. 4.1.1, the unnumbered
#' equation after eq 32), `2 b_MG - (b_MG^a + b_MG^b) / 2`, `a` and `b` the
#' estimates over periods `1..[T/2]` and `[T/2]+1..T`; each half is taken
#' as a panel of its own, so that the lags of `y` and of the averages of
#' its first rows come from inside the half, as "computed from the first
#' half of the available time period" reads, with the same `p`. `"rma"` is
#' recursive mean adjustment (sec. 4.1.2 and n. 13): `y` and each column
#' of `x` less the mean of the unit's values in the periods before,
#' `y_it - (t - 1)^-1 sum_{s<t} y_is`, a unit's first period dropped, and the
#' estimator applied to these series, so that the lags of `y` are lags of
#' the adjusted `y` and the averages are those of the adjusted values;
#' observed factors are not adjusted (Ditzen 2018, p. 589, leaves the
#' constant out). The paper states no variance for the corrected
#' estimators: `vmg` is the nonparametric mean group variance (eq 32)
#' applied to the corrected unit estimates (`bi`: `2 b_i - (b_i^a +
#' b_i^b) / 2` for the jackknife, the estimates on the adjusted series for
#' `"rma"`), and the uncorrected estimate and its variance are returned in
#' `mg0` and `vmg0` for comparison.
#'
#' **Five lessons.** Juodis and Reese (2026) draw five lessons for applied
#' work. (1) The averages are factor estimates measured with error, which
#' can bias the estimates and change their variance: results should be
#' checked against the number and kind of averages (`avg` keeps a subset),
#' and with slopes common to all units and `T` not small against `N` the
#' intervals of CCEP are not centred without a bias correction, which this
#' function does not make (sec. 3.1). `nfac` reports the number of factors
#' that the averages kept identify, by the two regularized procedures of
#' their appendix A.1 on a balanced panel: `er`, the eigenvalue ratio of
#' Ahn and Horenstein (2013) on the averages normalised by `Sigma^-1/2`,
#' `Sigma = sum_i (Z_i - Zbar)'(Z_i - Zbar) / (N T)` (Juodis 2022; steps
#' 1-3, `r_max` the averages less one); `thr`, the number of eigenvalues
#' of `Zbar'Zbar / T` at least 0.01 times the largest (Vogt et al. 2025;
#' steps 1-2); and `avg`, the averages kept; `NA` where undefined. More
#' averages than factors is the case in which the paper warns of extra
#' bias; the regularized estimators themselves are not made. (2) A two-way
#' fixed-effects estimator ([fepnl()]) may still be consistent and should
#' be compared and tested. `hb` is their test of heterogeneity bias in
#' CCEP (sec. 3.2, case study 2): `delta_i = (X_i'M_i X_i - N^-1 sum_j
#' X_j'M_j X_j) (b_i - b_MG)`, `HB = delta' V^-1 delta` with `delta` their
#' mean and `V = N^-2 sum delta_i delta_i'`, chi-squared on `k` degrees of
#' freedom (`stat`, `df`, `p`); a rejection says that CCEP is biased by
#' slopes that vary with the regressors. (3) Dynamic models need the
#' dynamics of every observed variable specified, and a bias correction is
#' recommended in short panels (sec. 3.3; `ylag`, `lags`, `bias` above).
#' (4) The CIPS panel unit-root test does not test what CCE estimation
#' needs; they recommend the CSB test of Pesaran et al. (2013) (sec. 3.4),
#' which this package does not provide. (5) There is no general recipe for
#' nonlinear models, where the average of `y` should not be among the
#' proxies (sec. 3.5); this function is linear. For cross-sectional
#' dependence in the residuals see [pcdtst()].
#'
#' **Fixed-T variance.** Brown, Schmidt and Wooldridge (2026, eqs 12-14
#' and Theorem 1) show that CCEP is the pooled least-squares estimate of
#' `y` on the regressors with the columns of `H` partialled out unit by
#' unit, and that time-varying regressors `d_t` added to it get
#' coefficients of exactly 0. Their fixed-T variance of CCEP, which counts
#' the sampling error of the averages (sec. 5.4, eqs 36-40, with `Psi = H`
#' and `q_i = vec(Z_i - Zbar)`, 0 in the columns of the observed factors),
#' is `pvar = "bsw"`: `A^-1 B A^-1 / N`, `A = N^-1 sum X_i'M X_i`, `B =
#' N^-1 sum s_i s_i'`, `s_i = X_i'M u_i - G'(I + K)[Psi (Psi'Psi)^-1 (x)
#' M] q_i`, `u_i = y_i - X_i b`, `G = N^-1 sum u_i (x) (X_i - Xbar)` and
#' `K` the commutation matrix (eq 33); it is for the static pooled
#' estimator with equal weights on a balanced panel, and refused
#' otherwise.
#'
#' @inheritParams gmmpnl
#' @param x The regressors: a numeric matrix or data frame, one row per row
#'   of `y`. For `plot()`, a `ccepnl` result.
#' @param fac `NULL` or observed common factors: a numeric matrix or data
#'   frame, one row per row of `y`, constant within each period.
#' @param cons `TRUE` (default) for a constant among the observed factors,
#'   unit-specific intercepts.
#' @param w `NULL` (equal weights `1 / N`, as Pesaran 2006, sec. 5.2, says
#'   is typical) or a positive weight of each row, constant within each
#'   unit; they are scaled to sum to 1 over the units.
#' @param avg `NULL` (default: the averages of `y` and of every column of
#'   `x`, Pesaran's choice) or the averages kept among the columns of `H`:
#'   `"y"` and names of columns of `x`.
#' @param ylag The lags of `y` among the regressors: a whole number from 0
#'   (default 0, the static model of Pesaran 2006).
#' @param lags `NULL` (default: 0 for the static model, the integer part of
#'   `T^(1/3)` with `ylag > 0`) or the lags of the averages: a whole number
#'   from 0.
#' @param bias `"none"` (default), `"jack"` (half-panel jackknife) or
#'   `"rma"` (recursive mean adjustment): the correction of the mean group
#'   estimate, as above.
#' @param pool `NULL` (default: CCEP for the static model only), `TRUE` or
#'   `FALSE`: whether to make CCEP; `TRUE` is refused with `ylag > 0` or
#'   `bias`.
#' @param pvar `"pesaran"` (default, eqs 67-69 of Pesaran 2006) or `"bsw"`
#'   (Brown, Schmidt and Wooldridge 2026): the variance of CCEP.
#' @return An object of class `ccepnl`, a list: `mg` and `pooled`, the
#'   CCEMG and CCEP estimates (`term`, `est`, `se`, `z`, `p`, `lo`, `hi`;
#'   `pooled` `NULL` when not made); `mg0`, the uncorrected mean group
#'   estimate with `bias` (`NULL` otherwise); `bi`, the unit estimates
#'   (units by regressors, corrected with `bias`); `vmg`, `vp` and `vmg0`,
#'   the variances; `hb`, the heterogeneity-bias test of CCEP (`NULL`
#'   without CCEP); `nfac`, the number of factors (`er`, `thr`, `avg`);
#'   `fit`, one row (`units`, `periods`, `k`, `h` the columns of `H`, `n`
#'   the rows of the regressions, `lags`, `ylag`); `w`, the weights by unit;
#'   `bias`; `pvar`; and `level`.
#' @references
#'   Pesaran, M. H. (2006). Estimation and inference in large heterogeneous
#'   panels with a multifactor error structure. *Econometrica*, 74(4),
#'   967-1012. \doi{10.1111/j.1468-0262.2006.00692.x}
#'
#'   Ditzen, J. (2018). Estimating dynamic common-correlated effects in
#'   Stata. *The Stata Journal*, 18(3), 585-617.
#'   \doi{10.1177/1536867X1801800306}
#'
#'   Su, L., Wang, F. and Wang, Y. (2026). Estimation and inference for
#'   unbalanced panel data models with interactive fixed effects. *Journal
#'   of Econometrics*, 255, 106222. \doi{10.1016/j.jeconom.2026.106222}
#'
#'   Juodis, A. and Reese, S. (2026). Five lessons for applied researchers
#'   from twenty years of common correlated effects estimation. *Journal of
#'   Econometrics*, 253, 106120. \doi{10.1016/j.jeconom.2025.106120}
#'
#'   Chudik, A. and Pesaran, M. H. (2015). Common correlated effects
#'   estimation of heterogeneous dynamic panel data models with weakly
#'   exogenous regressors. *Journal of Econometrics*, 188(2), 393-420.
#'   \doi{10.1016/j.jeconom.2015.03.007}
#'
#'   Brown, N. L., Schmidt, P. and Wooldridge, J. M. (2026). A simple
#'   reformulation of the common correlated effects model. *Econometric
#'   Reviews*, 45(5), 760-772. \doi{10.1080/07474938.2026.2615434}
#'
#'   Czarnowske, D. and Stammann, A. (2026). Inference in unbalanced panel
#'   data models with interactive fixed effects. *Econometric Reviews*,
#'   45(9), 1306-1337. \doi{10.1080/07474938.2026.2681519}
#' @examples
#' # 30 provinces over 12 years, one unobserved factor in y and x
#' set.seed(9)
#' d <- expand.grid(year = 1:12, prov = paste0("p", 1:30))
#' f <- rnorm(12)[d$year]
#' g <- rnorm(30, 1)[as.integer(substring(d$prov, 2))]
#' d$x <- g * f + rnorm(360)
#' d$y <- 1 + 0.7 * d$x + g * f + rnorm(360)
#' r <- ccepnl(d$y, d["x"], d$prov, d$year)
#' r
#' r$nfac
#' # a dynamic panel: y depends on its lag; jackknife-corrected mean group
#' set.seed(10)
#' tt <- 30
#' d <- expand.grid(year = 1:tt, prov = paste0("p", 1:40))
#' f <- rnorm(tt)[d$year]
#' g <- rnorm(40, 1)[as.integer(substring(d$prov, 2))]
#' d$x <- g * f + rnorm(nrow(d))
#' d$y <- 0
#' for (i in which(d$year > 1)) {
#'   d$y[i] <- 0.5 * d$y[i - 1] + 0.7 * d$x[i] + g[i] * f[i] + rnorm(1)
#' }
#' ccepnl(d$y, d["x"], d$prov, d$year, ylag = 1, bias = "jack")
#' @export
ccepnl <- function(y, x, unit, time, fac = NULL, cons = TRUE, w = NULL,
                   avg = NULL, ylag = 0, lags = NULL, bias = "none",
                   pool = NULL, pvar = "pesaran", level = 0.95) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(fac)) fac <- as.matrix(fac)
  .Call(C_clab_r_ccepnl, y, x, unit, time, fac, cons, w, avg, ylag, lags,
        bias, pool, pvar, level)
}

#' @export
print.ccepnl <- function(x, ...) {
  f <- x[["fit"]]
  b <- x[["bias"]]
  cat("ccepnl: common correlated effects, ", f[["units"]], " units, ",
      f[["periods"]], " periods",
      if (f[["ylag"]] > 0 || f[["lags"]] > 0) {
        paste0(", ", f[["ylag"]], " lags of y, ", f[["lags"]],
               " of the averages")
      },
      "\nmean group (CCEMG", if (!identical(b, "none")) paste0(", ", b),
      "):\n", sep = "")
  print(x[["mg"]], row.names = FALSE)
  if (!is.null(x[["pooled"]])) {
    cat("pooled (CCEP):\n")
    print(x[["pooled"]], row.names = FALSE)
    h <- x[["hb"]]
    cat("heterogeneity bias of CCEP: chi2(", h[["df"]], ") = ",
        format(h[["stat"]], digits = 4), ", p ", format(h[["p"]], digits = 3),
        "\n", sep = "")
  }
  invisible(x)
}

#' @rdname ccepnl
#' @inheritParams plotlrnz
#' @details The plot draws the CCEMG and, when made, CCEP estimates of each
#'   regressor with their intervals and returns what it drew.
#' @export
plot.ccepnl <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "ccepnl", "mg")
  a <- x$mg
  b <- x$pooled
  if (is.null(b)) b <- a[0, , drop = FALSE]
  .pdots(c(a$term, b$term), c(a$est, b$est), c(a$lo, b$lo), c(a$hi, b$hi),
         .pmain(main, "Common correlated effects"), "Coefficient", col, pal,
         ref = 0, grp = rep(c("CCEMG", "CCEP"), c(nrow(a), nrow(b))), ...)
}

# One line of a joint test of gmmpnl(): F or Wald chi-squared, with a tail.
.gwline <- function(w, tail) {
  paste0(if (w$test == "F") {
    paste0("F(", w$df1, ", ", w$df2, ")")
  } else {
    paste0("Wald chi2(", w$df1, ")")
  }, " = ", format(w$stat, digits = 5), ", p ", format(w$p, digits = 4),
  tail, "\n")
}
