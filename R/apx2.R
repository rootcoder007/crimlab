# Exposure and crime, the designs of 2019-2026: smoke instruments, runs of
# exposed days, heterogeneity, the case-crossover design, joint exposures,
# dose-response bins and effects by offence type. Each builds the columns
# of its design and fits it with ppmlfe(), ivfe() or ppmlcf() (R/apx.R).

#' Smoke instruments from plume data
#'
#' Builds the wildfire-smoke instruments of the recent papers from a daily
#' smoke-plume record by place: the presence of a plume; the presence
#' adjusted to the plumes that reached the ground; plume density classes;
#' smoke-attributable fine particles; and counts by month.
#'
#' The choices and their sources:
#'
#' * `adj = "raw"`: the plume indicator as given. Kircheis (2026, p. 4)
#'   instruments daily PM2.5 by "a daily binary indicator of predicted
#'   ground-level wildfire smoke presence" from a machine-learning product;
#'   Lee and Lee (2026, p. 6) use the same smoke days and, as a robustness
#'   check, NOAA's Hazard Mapping System plumes alone.
#' * `adj = "own"`: a plume day counts only when the place's PM2.5 is more
#'   than `z` standard deviations above its background, the mean of its
#'   days without a plume in the same season (Burkhardt et al. 2019,
#'   footnote 13: the own-county adjusted variable, 1.64 standard
#'   deviations, the 95th percentile under normality).
#' * `adj = "nbr"`: the same test made with the plume and the PM2.5 of the
#'   nearest neighbouring place (`nbr`), Burkhardt et al.'s preferred
#'   neighbour-adjusted variable (p. 11-12), which avoids using the place's
#'   own PM2.5 in its instrument.
#' * `seas`: the season of the background. Burkhardt et al. say "three
#'   month seasonal mean" without naming the months; the default is the
#'   meteorological seasons (December to February, ...) from a date column,
#'   pooled over years; give your own column to choose otherwise.
#' * `dens`: density classes (HMS "light", "medium", "heavy") become one
#'   indicator each. Lee and Lee (2026, footnote 4) warn that HMS density
#'   before 2018 has large measurement error, and use presence only.
#' * `spm`: smoke-attributable PM2.5, PM2.5 less the place-season
#'   background (`cf = "mean"`, Burkhardt et al.'s background, or
#'   `"median"`), on adjusted smoke days and at least 0; 0 on other days.
#'   This counterfactual is crimlab's construction from Burkhardt et al.'s
#'   background. `cf = "childs"` is instead the definition of Childs et
#'   al. (2022, equations 1-3, p. 13609): the background is the median
#'   PM2.5 of the same place and calendar month over the year and the
#'   years either side (three years) on days without smoke (`plume` 0),
#'   and smoke PM2.5 is PM2.5 above it on smoke days (`plume` 1), at least
#'   0, and 0 on other days; a place-month-year with no smoke-free days in
#'   its window gets `NA`, with a warning. It needs `day` as dates and
#'   ignores `seas` and `adj` (their smoke days are the plume column as
#'   given). Their machine-learned gridded smoke PM2.5 is a data product
#'   to be supplied as a column, not rebuilt here.
#' * `by`: with a group (place and month), the counts of smoke days and
#'   the means of the rest per group, Burkhardt et al.'s county-month
#'   instrument, "the number of days a county is affected by a fire smoke
#'   plume in a given month" (p. 10).
#'
#' Bao and Da (2026, sec. 3.3) keep only "Medium and Heavy plumes"
#' because "including Light plumes badly overstates surface-smoke days"
#' (give the density classes in `dens` and use their indicators), and
#' instead of a plume overhead count the days a county lies downwind of a
#' plume ([dwnexp()]).
#'
#' The instrument must move the exposure (a first stage, see [ivfe()]) and
#' change the outcome only through it: smoke may also keep people indoors
#' or darken the day (Lee and Lee 2026 discuss both), which the exclusion
#' restriction rules out by assumption.
#'
#' @param data A data frame, one row per place and day.
#' @param plume The name of the plume column (0 or 1).
#' @param unit,day The names of the place and day columns.
#' @param pm `NULL` or the name of the PM2.5 column (needed for `"own"`,
#'   `"nbr"` and `spm`).
#' @param adj `"raw"`, `"own"` or `"nbr"`.
#' @param nbr For `"nbr"`, the name of the column naming each row's nearest
#'   neighbouring place (a value of `unit`).
#' @param z The number of standard deviations above the background.
#' @param seas `NULL` (meteorological seasons of a date `day`) or the name
#'   of a season column.
#' @param dens `NULL` or the name of a density-class column (`NA` or `""`
#'   for no plume).
#' @param cf `"mean"`, `"median"` or `"childs"`: the background of `spm`.
#' @param by `NULL` or the name of a grouping column (place-month) to
#'   aggregate to.
#' @return A data frame, one row per row of `data` (or per `by` group):
#'   `smoke` (0 or 1, or the count of smoke days per group), `bg` (the
#'   background, with `pm`), `spm` (with `pm`), and one indicator (or count)
#'   per density class, `dens:<class>`.
#' @references
#'   Childs, M. L., Li, J., Wen, J., Heft-Neal, S., Driscoll, A., Wang, S.,
#'   Gould, C. F., Qiu, M., Burney, J. and Burke, M. (2022). Daily
#'   local-level estimates of ambient wildfire smoke PM2.5 for the
#'   contiguous US. *Environmental Science and Technology*, 56(19),
#'   13607-13621. \doi{10.1021/acs.est.2c02934}
#'
#'   Burkhardt, J., Bayham, J., Wilson, A., Berman, J. D., O'Dell, K.,
#'   Ford, B., Fischer, E. V. and Pierce, J. R. (2019). The relationship
#'   between monthly air pollution and violent crime across the United
#'   States. *Journal of Environmental Economics and Policy*, 9(2),
#'   188-205. \doi{10.1080/21606544.2019.1630014}
#'
#'   Kircheis, L. (2026). Wildfire smoke increases assaults: evidence from
#'   Seattle. *Environmental Research Letters*, 21, 044020.
#'   \doi{10.1088/1748-9326/ae436c}
#'
#'   Lee, G. and Lee, S. (2026). The impact of wildfire smoke exposure on
#'   crime. *Environmental and Resource Economics*, 89, 5.
#'   \doi{10.1007/s10640-025-01053-2}
#'
#'   Bao, J. and Da, Y. (2026). Wildfire smoke and crime: evidence from
#'   California. Working paper, Renmin University of China.
#' @seealso [ppmlcf()], [exprun()], [dwnexp()]
#' @examples
#' # simulated plumes and PM2.5 in two places over 120 days
#' set.seed(6)
#' d <- data.frame(place = rep(c("a", "b"), each = 120),
#'                 day = rep(as.Date("2024-06-01") + 0:119, 2))
#' d$plume <- rbinom(240, 1, 0.2)
#' d$pm <- rnorm(240, 8, 2) + d$plume * rbinom(240, 1, 0.5) * 15
#' d$near <- ifelse(d$place == "a", "b", "a")
#' s <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "own")
#' table(raw = d$plume, adjusted = s$smoke)
#' head(smkiv(d, "plume", "place", "day", pm = "pm", adj = "nbr", nbr = "near"))
#' @export
smkiv <- function(
  data,
  plume,
  unit,
  day,
  pm = NULL,
  adj = c("raw", "own", "nbr"),
  nbr = NULL,
  z = 1.64,
  seas = NULL,
  dens = NULL,
  cf = c("mean", "median", "childs"),
  by = NULL
) {
  f <- "smkiv"
  adj <- .ckmarg("adj", f)
  cf <- .ckmarg("cf", f)
  d <- .axdat(data, f)
  ck <- .axsmck(d, plume, unit, day, pm, seas, by, dens, adj, z, f)
  pl <- ck$pl
  z <- ck$z
  sv <- .axsmsea(d, day, seas, pm, f)
  key <- paste(as.character(d[[unit]]), sv, sep = "\r")
  out <- data.frame(smoke = pl)
  if (!is.null(pm)) {
    out <- .axsmpm(d, out, pl, key, pm, adj, nbr, z, cf, unit, day, f)
  }
  out <- .axsmdn(out, d, dens)
  if (!is.null(by)) return(.axsmby(out, d, by))
  out
}

#' Runs of exposed days
#'
#' Counts, for each place and day, how many days in a row the place has
#' been exposed (smoke present, pollution above a threshold), and makes the
#' indicators of the first, second, ... day of a run. `runfx()` fits a
#' model with them: the effect by day in the run, as Lee and Lee (2026,
#' equation 3, p. 12) estimate it, `Crime = b1 Day1 + b2 Day2 + b3 Day3+ +
#' ...`, which shows whether the effect fades as people adapt (it rises on
#' the first and second smoke days and falls from the third in their 21
#' cities).
#'
#' A run continues only to the next calendar day (`step`): a missing day
#' or an unexposed day ends it. `cap` sets the last indicator, which takes
#' that day and all later ones (3 in the paper: `run3p`). Lee and Lee use
#' Poisson pseudo-likelihood with ZIP-by-year, year-by-month, city-by-month
#' and weekend effects and errors clustered by city and day of the year
#' (p. 11 and table 3); [ppmlfe()] absorbs one or two sets of effects, so
#' give others as columns of `ctrl`.
#'
#' @inheritParams smkiv
#' @inheritParams dosebin
#' @param unit,day The names of the place and day columns (days numbers or
#'   dates).
#' @param thr A day is exposed when `var` is above `thr`.
#' @param cap The day of a run from which one indicator takes all, at most
#'   67108864 (2^26).
#' @param step The length of a day in the units of `day`.
#' @param est The estimator: [ppmlfe()] (the default) or [ivfe()] as least
#'   squares.
#' @param x An object of class `runfx`.
#' @return `exprun()`: a data frame with `run` (0 when not exposed) and the
#'   indicators `run1`, ..., `run<cap>p`. `runfx()`: a list of class
#'   `runfx`, `coef` (the rows of the indicators) and `fit`.
#' @references
#'   Lee, G. and Lee, S. (2026). The impact of wildfire smoke exposure on
#'   crime. *Environmental and Resource Economics*, 89, 5.
#'   \doi{10.1007/s10640-025-01053-2}
#' @seealso [smkiv()], [ppmlfe()]
#' @examples
#' d <- data.frame(area = 1, day = 1:8, smoke = c(0, 1, 1, 1, 1, 0, 1, 0))
#' cbind(d, exprun(d, "smoke", "area", "day"))
#' # simulated: crime up on the first two smoke days
#' set.seed(7)
#' d <- expand.grid(day = 1:120, area = 1:12)
#' d$smoke <- rbinom(nrow(d), 1, 0.3)
#' r <- exprun(d, "smoke", "area", "day")
#' d$crimes <- rpois(nrow(d), exp(1 + 0.2 * r$run1 + 0.1 * r$run2))
#' f <- runfx(d, "crimes", "smoke", "area", "day", fe = c("area", "day"))
#' f
#' plot(f)
#' @export
exprun <- function(data, var, unit, day, thr = 0, cap = 3, step = 1) {
  f <- "exprun"
  d <- .axdat(data, f)
  .axcols(d, c(var, unit, day), "var, unit, day", f, n = 3L)
  thr <- .cknum(thr, "thr", f)
  cap <- .cknum(cap, "cap", f, min = 1, whole = TRUE)
  if (cap > 2^26) .ckstop(f, "`cap` must be at most 67108864")
  step <- .cknum(step, "step", f, positive = TRUE)
  tv <- d[[day]]
  if (!is.numeric(tv) && !inherits(tv, "Date")) {
    .ckstop(f, "column `%s` must hold numbers or dates", day)
  }
  e <- as.numeric(.axmat(d, var, f)[, 1L] > thr)
  run <- .Call(C_clab_r_runlen, .axcode(d[[unit]])$c, as.numeric(tv), e, step)
  out <- data.frame(run = run)
  for (k in seq_len(cap)) {
    nm <- if (k < cap) paste0("run", k) else paste0("run", k, "p")
    out[[nm]] <- as.numeric(if (k < cap) run == k else run >= k)
  }
  out
}

#' @rdname exprun
#' @export
runfx <- function(data, y, var, unit, day, thr = 0, cap = 3, step = 1,
                  ctrl = NULL, est = ppmlfe, ...) {
  f <- "runfx"
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  r <- .clas(f, exprun(data, var, unit, day, thr = thr, cap = cap,
                       step = step))
  rc <- names(r)[-1L]
  d <- data
  for (k in rc) d[[k]] <- r[[k]]
  fit <- .clas(f, est(data = d, y = y, x = c(rc, ctrl), ...))
  structure(
    list(coef = fit$coef[fit$coef$term %in% rc, , drop = FALSE], fit = fit),
    class = "runfx"
  )
}

#' @rdname exprun
#' @export
print.runfx <- function(x, ...) {
  cat("runfx: the effect by day of a run of exposed days\n")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' Effects by a moderator, with the joint test of equality
#'
#' Estimates the effect of an exposure separately for each level of a
#' moderator (household wealth, income quantile, rural or urban, offence
#' context) and tests whether the effects are equal.
#'
#' The choices and their sources: `how = "split"` fits the estimator on
#' each level's rows alone, as Kishore, Tsaneva and Balakrishnan (2026,
#' section 4.7 and table 5: poor and non-poor households, rural and urban)
#' and Bondy, Roth and Sager (2020, figure 5: house-price quantiles) do;
#' every coefficient and fixed effect then differs by level. The test of
#' equality treats the levels' estimates as independent, which holds when
#' no cluster spans two levels. `how = "inter"` (the default) fits one
#' model with the exposure interacted with each level's indicator (and, for
#' an instrumental-variable estimator, each instrument too), the other
#' coefficients and the fixed effects common; Bondy, Roth and Sager report
#' that this "joint regression including quantile interactions" gives
#' similar estimates (footnote 29). Its test uses the fit's covariance, so
#' it allows for clusters spanning levels. Neither paper reports a test of
#' equality; it is the Wald test of `b_1 = b_2 = ... = b_K` on `K - 1`
#' degrees of freedom.
#'
#' @inheritParams dosebin
#' @param var The name of the exposure (a regressor `x` or endogenous
#'   `endo` of the estimator).
#' @param mod The name of the moderator column (labels).
#' @param est The estimator: [ppmlfe()], [ivfe()] or [ppmlcf()].
#' @param how `"inter"` or `"split"`.
#' @param ... The other arguments of `est`. For `plot()`, passed to the
#'   base-graphics call.
#' @param x An object of class `hetfx`.
#' @return A list of class `hetfx`: `tab`, one row per level (`level`,
#'   `est`, `se`, `lo`, `hi`, `p`, `n`); `test`, the Wald test of equality
#'   (`stat`, `df`, `p`); `how`.
#' @references
#'   Kishore, S., Tsaneva, M. and Balakrishnan, U. (2026). Air pollution and
#'   intimate partner violence in India. *Review of Economics of the
#'   Household*. \doi{10.1007/s11150-026-09851-3}
#'
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#' @seealso [typefx()], [ppmlfe()]
#' @examples
#' # simulated: pollution matters in poor areas only
#' set.seed(8)
#' d <- expand.grid(day = 1:60, area = 1:20)
#' d$poor <- ifelse(d$area <= 10, "poor", "not poor")
#' d$pm <- rnorm(nrow(d), 20, 5)
#' d$crimes <- rpois(nrow(d), exp(0.5 + 0.04 * d$pm * (d$poor == "poor")))
#' h <- hetfx(d, "pm", "poor", y = "crimes", x = "pm", fe = c("area", "day"))
#' h
#' plot(h)
#' @export
hetfx <- function(data, var, mod, est = ppmlfe, how = c("inter", "split"),
                  ...) {
  f <- "hetfx"
  how <- .ckmarg("how", f)
  d <- .axdat(data, f)
  .axcols(d, c(var, mod), "var, mod", f, n = 2L)
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  mv <- as.character(d[[mod]])
  lv <- sort(unique(mv))
  if (length(lv) < 2L) {
    .ckstop(f, "column `%s` has one level: nothing to compare", mod)
  }
  a <- list(...)
  h <- if (how == "split") {
    .clas(f, .axhtspl(d, est, a, var, mv, lv, f))
  } else {
    .clas(f, .axhtint(d, est, a, var, mv, lv))
  }
  tab <- h$tab
  b <- h$b
  v <- h$v
  k <- length(b)
  rr <- cbind(diag(k - 1L), 0) - cbind(0, diag(k - 1L))
  structure(
    list(tab = tab, test = .axwald(b, v, rr), how = how),
    class = "hetfx"
  )
}

#' @rdname hetfx
#' @export
print.hetfx <- function(x, ...) {
  cat(
    "hetfx: the effect by level (",
    if (x[["how"]] == "inter") "interactions" else "separate fits",
    ")\n",
    sep = ""
  )
  print(x[["tab"]], row.names = FALSE, digits = 4)
  t <- x[["test"]]
  cat("Wald test of equal effects: ", format(t$stat, digits = 4), " on ", t$df,
      " df, p = ", format(t$p, digits = 3), "\n", sep = "")
  invisible(x)
}

#' The time-stratified case-crossover design
#'
#' Compares each day's count with the counts of the same place on the
#' other days of its stratum, the case-crossover design for counts: a
#' Poisson model with one fixed effect per stratum, which gives the same
#' coefficients as the conditional Poisson model. Ju et al. (2026,
#' equations 1 and 2, p. 7) use the space-time stratum "one region nested
#' in the same weekday of same month from the same year" (`strat =
#' "ymw"`, the default), a control function for wildfire PM2.5 driven by
#' wind speed, the 2-day moving average of exposure (`ma`), and errors
#' clustered by calendar year-month (the default `clus` here).
#'
#' Armstrong, Gasparrini and Tobias (2014, p. 3-4 and table 4) show that
#' the Poisson model with stratum indicators gives the same estimates as
#' the conditional Poisson and the conditional logistic case-crossover
#' model, and that the Poisson forms, unlike the conditional logistic one,
#' allow for overdispersion (quasi-Poisson, the dispersion estimated from
#' Pearson's chi-squared: `se = "quasi"` of [ppmlfe()], passed through
#' `...`), for autocorrelation and for varying denominators (an offset:
#' `expo`). The clustered errors used by default here also allow for
#' overdispersion and for correlation within a year-month. Their
#' adjustment for autocorrelation (Brumback's method, in their additional
#' file) is not built: the paper does not state it.
#'
#' Han et al. (2026, Methods, p. 2) run the design on 483,718 crime events
#' in China, 2013-19: control days the same weekday of the same calendar
#' month (`strat = "ymw"`), conditional logistic regression (equivalent to
#' the Poisson with stratum effects on daily counts by place, as above),
#' an indicator of public holidays, natural cubic splines of temperature
#' (6 degrees of freedom) and relative humidity (3) (`ns`), single-day lags
#' 0 to 6 (`lagv`, `lag`; their main table is lag 0), exposure-response
#' curves by a natural spline of the pollutant with 3 degrees of freedom
#' (give the pollutant in `ns` and plot its fitted terms), effect
#' modification by stratified fits compared by two-sample z-tests (fit
#' each subset, or [hetfx()] with `how = "split"`), and two-pollutant
#' models (two exposures in `x`), reporting percent changes per
#' interquartile range, `100 (exp(b IQR) - 1)`. Their weighted quantile sum
#' regression and air quality score are not built.
#'
#' The choices: `strat = "ymw"` (place, year, month and weekday) or `"ym"`
#' (place, year and month, the plain time-stratified design); `ma`, the
#' columns to replace by their mean over the day and the `k - 1` days
#' before (2 in Ju et al. and in Heo et al. 2026, "lag01"; a day whose
#' earlier days are missing is dropped); `est`, [ppmlfe()] for an
#' exposure taken as exogenous or [ppmlcf()] for the instrumented design of
#' Ju et al. (give `endo` and `inst`).
#'
#' Deaths or offences one by one. Rahman et al. (2023, sec. 2) fit the
#' individual case-crossover design by conditional logistic regression: for
#' each death the day it happened against "the first, second, and fourth
#' Tuesdays" of the month when it fell on the third, the exposure that of
#' the decedent's census tract on each day, each lag from 0 to 7 in its own
#' model, relative humidity as a natural spline of 3 degrees of freedom.
#' When every decedent of a tract shares the tract's exposure, the
#' conditional logistic likelihood of the individuals is the conditional
#' Poisson likelihood of the tract's daily counts: the individuals of one
#' tract, month and weekday share their referent days and exposures, so
#' that design is this one with `unit` the tract, `strat = "ymw"`, `lagv`
#' and `lag` for the lag and `ns = c(rh = 3)`, its coefficients the log
#' odds ratios.
#'
#' @inheritParams smkiv
#' @inheritParams jntx
#' @param unit The name of the place column.
#' @param day The name of the day column (dates).
#' @param strat `"ymw"` or `"ym"`.
#' @param ma `NULL` or the names of columns to average over days.
#' @param k The days of the moving average.
#' @param lagv `NULL` or the names of columns to replace by their value
#'   `lag` days before (single-day lags).
#' @param lag The lag in days for `lagv`.
#' @param ns `NULL` or degrees of freedom named by column, for example
#'   `c(temp = 6, rh = 3)`: natural cubic splines of those columns join
#'   the controls `x` (named `"<col>:ns<j>"`).
#' @param est The estimator, [ppmlfe()] or [ppmlcf()].
#' @param ... The other arguments of `est` (`x`, or `endo` and `inst`;
#'   `clus` to change the clustering).
#' @return The fit of `est`, with the stratum column `.strat` and the
#'   year-month column `.ym` added to the data it used.
#' @references
#'   Ju, K., Xu, R., Huang, W., Zhang, Y., Ye, T., Yu, W., Yang, Z., Chen,
#'   G., Liu, Y., Yu, P., Zhou, S., Wu, Y., Wen, B., Xu, Z., Saldiva, P. H.
#'   N., Coelho, M. S. Z. S., Song, J., Li, S. and Guo, Y. (2026). Causal
#'   effects of wildfire PM2.5 on hospital costs and length of stay in
#'   Brazil. *Nature Communications*. \doi{10.1038/s41467-026-75157-8}
#'
#'   Heo, S., Berman, J. D. and Bell, M. L. (2026). Multi-city study of the
#'   joint effects of ozone and temperature exposure on violent crime during
#'   warm months in the United States. *Environmental Research*,
#'   310, 125969. \doi{10.1016/j.envres.2026.125969}
#'
#'   Armstrong, B. G., Gasparrini, A. and Tobias, A. (2014). Conditional
#'   Poisson models: a flexible alternative to conditional logistic case
#'   cross-over analysis. *BMC Medical Research Methodology*, 14, 122.
#'
#'   Burkhardt, J., Bayham, J., Wilson, A., Carter, E., Berman, J. D.,
#'   O'Dell, K., Ford, B., Fischer, E. V. and Pierce, J. R. (2019). The
#'   effect of pollution on crime: evidence from data on particulate matter
#'   and ozone. *Journal of Environmental Economics and Management*,
#'   98, 102267. \doi{10.1016/j.jeem.2019.102267}
#'
#'   Han, D., Guo, Y., Jiang, Y., Niu, Y., Chen, R., Kan, H. and others
#'   (2026). Air pollution and criminal behavior: a nationwide
#'   case-crossover analysis in China. *International Journal of
#'   Epidemiology*, 55(5), dyag214. \doi{10.1093/ije/dyag214}
#'
#'   Rahman, M. M., Lorenzo, M., Ban-Weiss, G., Hasan, Z., Azzouz, M.,
#'   Eckel, S. P., Conti, D. V., Lurmann, F., Schlaerth, H., Johnston, J.,
#'   Ko, J., Palinkas, L., Hurlburt, M., Silva, S., Gauderman, W. J.,
#'   McConnell, R. and Garcia, E. (2023). Ambient temperature and air
#'   pollution associations with suicide and homicide mortality in
#'   California: a statewide case-crossover study. *Science of the Total
#'   Environment*, 874, 162462. \doi{10.1016/j.scitotenv.2023.162462}
#' @seealso [ppmlfe()], [ppmlcf()], [jntx()]
#' @examples
#' set.seed(9)
#' d <- expand.grid(day = as.Date("2023-01-01") + 0:364, area = 1:4)
#' d$pm <- rnorm(nrow(d), 15, 4)
#' d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$pm))
#' ccross(d, "crimes", "area", "day", ma = "pm", x = "pm")
#'
#' # shaped as Burkhardt et al. (2019, JEEM, equation 1, p. 5): county-day
#' # counts on PM2.5 and ozone together, county-by-year-by-month-by-weekday
#' # effects, errors clustered by county (simulated)
#' d$o3 <- rnorm(nrow(d), 0.03, 0.01)
#' ccross(d, "crimes", "area", "day", x = c("pm", "o3"), clus = "area")
#'
#' # shaped as Han et al. (2026, Methods p. 2 and table 1): city-day
#' # counts, same-weekday-same-month strata, a holiday indicator, natural
#' # splines of temperature (6 df) and humidity (3 df), NO2 at lag 0, the
#' # percent change per interquartile range (simulated)
#' d$no2 <- rnorm(nrow(d), 30, 10)
#' d$temp <- 15 + 10 * sin(2 * pi * as.numeric(d$day) / 365) +
#'   rnorm(nrow(d))
#' d$rh <- runif(nrow(d), 30, 90)
#' d$hol <- as.numeric(format(d$day, "%m-%d") %in% c("01-01", "10-01"))
#' d$crimes <- rpois(nrow(d), exp(1 + 0.004 * d$no2))
#' h <- ccross(d, "crimes", "area", "day", x = c("no2", "hol"),
#'             ns = c(temp = 6, rh = 3), clus = "area")
#' b <- h$coef[h$coef$term == "no2", ]
#' 100 * (exp(c(b$est, b$lo, b$hi) * IQR(d$no2)) - 1)
#' plot(h, terms = "no2")
#' @export
ccross <- function(data, y, unit, day, strat = c("ymw", "ym"), ma = NULL, k = 2,
                   lagv = NULL, lag = 0, ns = NULL, est = ppmlfe, ...) {
  f <- "ccross"
  strat <- .ckmarg("strat", f)
  d <- .axdat(data, f)
  ck <- .axccck(d, y, unit, day, ma, est, k, lagv, lag, ns, f)
  d <- .clas(f, .axccma(d, unit, day, ma, ck$k, lagv, ck$lag))
  ym <- format(d[[day]], "%Y-%m")
  d$.strat <- paste(
    as.character(d[[unit]]),
    ym,
    if (strat == "ymw") as.POSIXlt(d[[day]])$wday else "",
    sep = "\r"
  )
  d$.ym <- ym
  a <- list(...)
  sp <- .axccns(d, a, ns, f)
  d <- sp$d
  a <- sp$a
  if (is.null(a[["clus"]])) a[["clus"]] <- ".ym"
  .clas(f, do.call(est, c(list(data = d, y = y, fe = ".strat"), a)))
}

#' Two exposures together: their interaction and the surface of effects
#'
#' Estimates how two exposures (heat and ozone, smoke and heat) act
#' together on a count, as a surface of rate ratios over both, the effect
#' of one at chosen levels of the other, and a test of the interaction.
#'
#' Heo, Berman and Bell (2026, p. 3) fit, for each city, a quasi-Poisson
#' model with a tensor-product smooth of the two-day moving averages
#' (lag01) of temperature and ozone, basis dimensions four and three, and
#' report the heat effect (the 90th against the 10th percentile of
#' temperature) at low, middle and high ozone (the 5th, 50th and 95th
#' percentiles), pooled over cities by a meta-analysis. The choices here:
#' `form = "spline"` (the default) uses the tensor product of two
#' restricted cubic spline bases (Harrell), `k1` and `k2` knots at fixed
#' quantiles (5 and 4 give the paper's dimensions four and three), without
#' the penalty of the paper's smooth, so the fit is that of a fixed basis,
#' not of 'mgcv'; `form = "linear"` uses `x1`, `x2` and `x1 x2`, the
#' common product-term model. The meta-analysis across cities is not
#' done here: fit each city and pool the `cond` rows. Berman et al. (2019)
#' estimate the risk of violent and property crime from daily PM2.5, ozone
#' and temperature in some 300 US counties with "a two-stage hierarchical
#' time series model", the counties' series pooled in the second stage,
#' the same pattern. Use [ccross()]'s
#' `ma` for the lag01 averages, and the fixed effects of [ppmlfe()] for
#' the paper's spline of calendar time and day of week (or give them in
#' `...`).
#'
#' @inheritParams dosebin
#' @param y The name of the count column.
#' @param x1,x2 The names of the two exposures.
#' @param form `"spline"` or `"linear"`.
#' @param k1,k2 The knots of the spline bases (3 to 7).
#' @param p1 The two quantiles of `x1` compared (`cond`).
#' @param p2 The quantiles of `x2` at which they are compared.
#' @param ref `NULL` (the medians) or the reference point of the surface.
#' @param ngrid The points of the surface along each exposure, at most
#'   8192 (`ngrid` times `ngrid` at most 2^26).
#' @param est The estimator, [ppmlfe()] by default.
#' @param x An object of class `jntx`.
#' @return A list of class `jntx`: `cond`, the log rate ratio (and `rr`)
#'   of `x1` at `p1[2]` against `p1[1]` at each quantile of `x2`, with
#'   `se`, `lo`, `hi`; `surf`, the grid (`x1`, `x2`, `est`, `se`, `rr`)
#'   relative to `ref`; `test`, the Wald test that every interaction term
#'   is 0; `fit`.
#' @references
#'   Heo, S., Berman, J. D. and Bell, M. L. (2026). Multi-city study of the
#'   joint effects of ozone and temperature exposure on violent crime during
#'   warm months in the United States. *Environmental Research*,
#'   310, 125969. \doi{10.1016/j.envres.2026.125969}
#'
#'   Berman, J., Burkhardt, J., Bayham, J., Carter, E. and Wilson, A.
#'   (2019). Acute air pollution exposure, temperature, and the risk of
#'   violent behavior in the United States. Conference abstract,
#'   *Environmental Epidemiology* (ISEE 2019 abstracts).
#'
#'   Harrell, F. E. (2015). *Regression Modeling Strategies*, 2nd ed.
#'   Springer. \doi{10.1007/978-3-319-19425-7}
#' @seealso [ccross()], [dosebin()], [ppmlfe()]
#' @examples
#' set.seed(10)
#' d <- expand.grid(day = 1:150, city = 1:10)
#' d$temp <- rnorm(nrow(d), 25, 4)
#' d$o3 <- 0.6 * d$temp + rnorm(nrow(d), 30, 6)
#' d$crimes <- rpois(
#'   nrow(d),
#'   exp(1 + 0.02 * d$temp + 0.001 * d$temp * d$o3 / 10)
#' )
#' j <- jntx(d, "crimes", "temp", "o3", fe = "city")
#' j
#' plot(j)
#' @export
jntx <- function(data, y, x1, x2, form = c("spline", "linear"), k1 = 5, k2 = 4,
                 p1 = c(
                   0.1,
                   0.9
                 ), p2 = c(0.05, 0.5, 0.95), ref = NULL, ngrid = 15,
                 ctrl = NULL, est = ppmlfe, ...) {
  f <- "jntx"
  form <- .ckmarg("form", f)
  d <- .axdat(data, f)
  ck <- .axjxck(d, y, x1, x2, est, k1, k2, p1, p2, ngrid, f)
  k1 <- ck$k1
  k2 <- ck$k2
  p1 <- ck$p1
  p2 <- ck$p2
  ngrid <- ck$ngrid
  if (!is.null(ref)) {
    ref <- .cknum(ref, "ref", f, scalar = FALSE)
    if (length(ref) != 2L) .ckstop(f, "`ref` must be two numbers, x1 and x2")
  }
  a <- .axmat(d, c(x1, x2), f)
  t1 <- .axknots(a[, 1L], k1)
  t2 <- .axknots(a[, 2L], k2)
  if (form == "spline" && (anyDuplicated(t1) || anyDuplicated(t2))) {
    .ckstop(f, "an exposure has too few distinct values for its knots")
  }
  basis <- function(u, v) .axjxbas(u, v, form, t1, t2)
  bm <- basis(a[, 1L], a[, 2L])
  nb1 <- if (form == "linear") 1L else k1 - 1L
  nb2 <- if (form == "linear") 1L else k2 - 1L
  nm <- paste0(".j", seq_len(ncol(bm)))
  for (j in seq_along(nm)) d[[nm[j]]] <- bm[, j]
  r <- .clas(f, est(data = d, y = y, x = c(nm, ctrl), ...))
  b <- r$coef$est[match(nm, r$coef$term)]
  v <- .axvcov(r)[nm, nm, drop = FALSE]
  if (is.null(ref)) {
    ref <- c(stats::median(a[, 1L]), stats::median(a[, 2L]))
  }
  rel <- function(u, w, u0, w0) {
    db <- basis(u, w) - basis(u0, w0)
    e <- as.numeric(db %*% b)
    s <- sqrt(pmax(rowSums((db %*% v) * db), 0))
    list(e = e, s = s)
  }
  o <- .axjxout(a, p1, p2, ref, ngrid, rel)
  ii <- (nb1 + nb2 + 1L):ncol(bm)
  rr <- diag(ncol(bm))[ii, , drop = FALSE]
  structure(
    list(cond = o$cond, surf = o$surf, test = .axwald(b, v, rr), fit = r,
         form = form),
    class = "jntx"
  )
}

#' @rdname jntx
#' @export
print.jntx <- function(x, ...) {
  cat(
    "jntx: the effect of x1 (",
    x[["form"]],
    " surface) at levels of x2\n",
    sep = ""
  )
  print(x[["cond"]], row.names = FALSE, digits = 4)
  t <- x[["test"]]
  cat("Wald test of no interaction: ", format(t$stat, digits = 4), " on ", t$df,
      " df, p = ", format(t$p, digits = 3), "\n", sep = "")
  invisible(x)
}

#' Exposure in bins: dose-response and tests of non-linearity
#'
#' Replaces a continuous exposure by indicators of its bins, the
#' reference bin left out, fits the model, and tests whether the bins
#' matter and whether their effects follow a straight line.
#'
#' Chin, Lurbe and Shakya (2026, section 5.2, p. 12) replace the log AQI by
#' indicators of its quintiles, the first left out (`nq = 5`, the
#' default), and count the days above public-health thresholds of 100 and
#' 150 (p. 13; give fixed breaks in `brk`, or aggregate days above a
#' threshold yourself); their main equation (2) is least squares of the log
#' count on the log exposure with police-agency-by-month-year effects and
#' errors clustered by district and season. They report no test of
#' non-linearity. The tests here: `all`, the Wald test that every bin's
#' coefficient is 0; `linear`, the Wald test that the coefficients are
#' proportional to the distance of each bin's median exposure from the
#' reference bin's (a straight line through the reference), on the number
#' of bins less two degrees of freedom.
#'
#' @param data A data frame.
#' @param y The name of the outcome column.
#' @param var The name of the exposure column.
#' @param nq The number of quantile bins (when `brk` is `NULL`), at most
#'   67108864 (2^26).
#' @param brk `NULL` or increasing breaks (right-closed; the lowest
#'   included).
#' @param ref The reference bin, counted from 1.
#' @param ctrl `NULL` or the names of other regressors.
#' @param est The estimator, [ppmlfe()] or [ivfe()] as least squares.
#' @param ... The other arguments of `est` (`fe`, `clus`, `se`, ...).
#'   For `plot()`, passed to the base-graphics call.
#' @param x An object of class `dosebin`.
#' @return A list of class `dosebin`: `tab`, one row per bin (`bin`, `mid`,
#'   the median exposure, `n`, `est`, `se`, `lo`, `hi`, `p`; the reference
#'   with 0); `test`, the two Wald tests; `fit`.
#' @references
#'   Chin, S., Lurbe, S. and Shakya, P. (2026). Air pollution and crime:
#'   heterogeneity by offense types and exposure intensity. SSRN working
#'   paper. \doi{10.2139/ssrn.6906621}
#' @seealso [jntx()], [typefx()]
#' @examples
#' set.seed(11)
#' d <- expand.grid(day = 1:100, area = 1:10)
#' d$aqi <- rgamma(nrow(d), 4, 0.1)
#' d$crimes <- rpois(nrow(d), exp(1 + 0.3 * (d$aqi > 50)))
#' b <- dosebin(d, "crimes", "aqi", fe = c("area", "day"))
#' b
#' plot(b)
#' @export
dosebin <- function(data, y, var, nq = 5, brk = NULL, ref = 1, ctrl = NULL,
                    est = ppmlfe, ...) {
  f <- "dosebin"
  d <- .axdat(data, f)
  .axcols(d, c(y, var), "y, var", f, n = 2L)
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  xv <- .axmat(d, var, f)[, 1L]
  bb <- .axdbbin(xv, nq, brk, var, f)
  bn <- bb$bn
  lv <- bb$lv
  nbin <- bb$nbin
  ref <- .cknum(ref, "ref", f, min = 1, max = nbin, whole = TRUE)
  nm <- paste0("bin:", lv)
  use <- seq_len(nbin)[-ref]
  for (j in use) d[[nm[j]]] <- as.numeric(bn == lv[j])
  r <- .clas(f, est(data = d, y = y, x = c(nm[use], ctrl), ...))
  cf <- r$coef[match(nm[use], r$coef$term), , drop = FALSE]
  mid <- as.numeric(tapply(xv, bn, stats::median))
  tab <- data.frame(
    bin = lv,
    mid = mid,
    n = as.integer(table(bn)),
    est = 0,
    se = NA_real_,
    lo = NA_real_,
    hi = NA_real_,
    p = NA_real_
  )
  tab[
    use,
    c("est", "se", "lo", "hi", "p")
  ] <- cf[c("est", "se", "lo", "hi", "p")]
  b <- cf$est
  v <- .axvcov(r)[nm[use], nm[use], drop = FALSE]
  k <- length(b)
  tst <- .axdbtst(b, v, mid, use, ref, k)
  structure(list(tab = tab, test = tst[c("test", "stat", "df", "p")], fit = r),
            class = "dosebin")
}

#' @rdname dosebin
#' @export
print.dosebin <- function(x, ...) {
  cat("dosebin: the effect of each bin against the reference\n")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  print(x[["test"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' Effects by offence type in one call, with a multiple-testing correction
#'
#' Fits the same model to several outcomes (assault, robbery, burglary,
#' property crime, ...) and adjusts the p-values for the number of
#' outcomes tested.
#'
#' Bondy, Roth and Sager (2020, figures 6 and 7) and Chin, Lurbe and
#' Shakya (2026) estimate the effect for each crime type separately and
#' report unadjusted tests; with many types, some will look significant by
#' chance. `adj` chooses the correction, through [stats::p.adjust()]:
#' `"holm"` (the default; it controls the chance of any false rejection
#' and is uniformly more powerful than Bonferroni), `"bonferroni"`,
#' `"BH"` (Benjamini-Hochberg, the false discovery rate) or `"none"`.
#' Corrections that use the correlation of the estimates (Romano-Wolf) are
#' not offered.
#'
#' @inheritParams dosebin
#' @param ys The names of the outcome columns, one fit each.
#' @param term The coefficient reported (the exposure).
#' @param est The estimator.
#' @param adj `"holm"`, `"bonferroni"`, `"BH"` or `"none"`.
#' @param ... The other arguments of `est` (`x`, `fe`, ...).
#'   For `plot()`, passed to the base-graphics call.
#' @param x An object of class `typefx`.
#' @return A list of class `typefx`: `tab`, one row per outcome (`y`,
#'   `est`, `se`, `lo`, `hi`, `p`, `padj`, `n`); `adj`.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Chin, S., Lurbe, S. and Shakya, P. (2026). Air pollution and crime:
#'   heterogeneity by offense types and exposure intensity. SSRN working
#'   paper. \doi{10.2139/ssrn.6906621}
#'
#'   Holm, S. (1979). A simple sequentially rejective multiple test
#'   procedure. *Scandinavian Journal of Statistics*, 6(2), 65-70.
#' @seealso [hetfx()], [ppmlfe()]
#' @examples
#' set.seed(12)
#' d <- expand.grid(day = 1:60, area = 1:12)
#' d$pm <- rnorm(nrow(d), 20, 5)
#' d$assault <- rpois(nrow(d), exp(0.5 + 0.03 * d$pm))
#' d$theft <- rpois(nrow(d), 2)
#' r <- typefx(d, c("assault", "theft"), "pm", x = "pm",
#'             fe = c("area", "day"))
#' r
#' plot(r)
#' @export
typefx <- function(data, ys, term, est = ppmlfe,
                   adj = c("holm", "bonferroni", "BH", "none"), ...) {
  f <- "typefx"
  adj <- .ckmarg("adj", f)
  d <- .axdat(data, f)
  .axcols(d, ys, "ys", f)
  if (!is.character(term) || length(term) != 1L) {
    .ckstop(f, "`term` must name one coefficient")
  }
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  rows <- lapply(ys, function(v) {
    r <- .clas(f, est(data = d, y = v, ...))
    cf <- .axrow(r, term, f)
    data.frame(
      y = v,
      est = cf$est,
      se = cf$se,
      lo = cf$lo,
      hi = cf$hi,
      p = cf$p,
      n = r$n
    )
  })
  tab <- do.call(rbind, rows)
  tab$padj <- stats::p.adjust(tab$p, method = adj)
  structure(
    list(
      tab = tab[c("y", "est", "se", "lo", "hi", "p", "padj", "n")],
      adj = adj
    ),
    class = "typefx"
  )
}

#' @rdname typefx
#' @export
print.typefx <- function(x, ...) {
  cat(
    "typefx: one fit per outcome, p-values adjusted by ",
    x[["adj"]],
    "\n",
    sep = ""
  )
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}

# A Wald test of R b = 0: the statistic, its degrees of freedom and the
# chi-square p-value.
# An estimator: a function with a `data` argument, as ppmlfe() and ivfe()
# have (sum(), say, would be called with the rows as its numbers).
.axest <- function(est) {
  is.function(est) && "data" %in% names(formals(args(est)))
}

.axwald <- function(b, v, r) {
  rb <- r %*% b
  m <- r %*% v %*% t(r)
  # a cluster-robust covariance has rank below the number of clusters, so
  # it cannot test more restrictions than that: no statistic, not a wrong
  # one from a matrix that is singular up to rounding
  if (qr(m)$rank < nrow(r)) {
    warning(
      "the joint test of ", nrow(r), " restrictions is not computed: the ",
      "covariance has rank ", qr(m)$rank, " (too few clusters for it).",
      call. = FALSE
    )
    return(data.frame(stat = NA_real_, df = nrow(r), p = NA_real_))
  }
  w <- as.numeric(t(rb) %*% solve(m, rb))
  data.frame(
    stat = w,
    df = nrow(r),
    p = stats::pchisq(w, nrow(r), lower.tail = FALSE)
  )
}

# The covariance of a fit: vcov, or from the bootstrap draws.
.axvcov <- function(r) {
  if (!is.null(r[["vcov"]])) return(r[["vcov"]])
  stats::cov(r[["draws"]], use = "complete.obs")
}

# The row of term in a fit's coefficients, or an error.
.axrow <- function(r, term, fn) {
  cf <- r[["coef"]]
  if (!is.data.frame(cf) || !term %in% cf[["term"]]) {
    .ckstop(fn, "the estimator returned no coefficient for `%s`", term)
  }
  cf[cf[["term"]] == term, , drop = FALSE]
}

# Replace var by new names in the estimator's arguments x and endo, and
# each instrument by its interactions (inst).
.axrepl <- function(a, var, new, inst = NULL) {
  for (k in c("x", "endo")) {
    if (!is.null(a[[k]]) && var %in% a[[k]]) {
      a[[k]] <- c(setdiff(a[[k]], var), new)
    }
  }
  if (!is.null(inst)) a[["inst"]] <- inst
  a
}

# A restricted cubic spline basis (Harrell): x and k - 2 cubic terms, with
# knots t at quantiles of x.
.axrcs <- function(x, t) {
  k <- length(t)
  if (k < 3L) return(matrix(x))
  p <- function(u) pmax(u, 0)^3
  m <- cbind(x, vapply(seq_len(k - 2L), function(j) {
    (p(x - t[j]) - p(x - t[k - 1L]) * (t[k] - t[j]) / (t[k] - t[k - 1L]) +
       p(x - t[k]) * (t[k - 1L] - t[j]) / (t[k] - t[k - 1L])) / (t[k] - t[1L])^2
  }, numeric(length(x))))
  unname(m)
}

.axknots <- function(x, k) {
  pr <- switch(
    as.character(k),
    `3` = c(0.1, 0.5, 0.9),
    `4` = c(0.05, 0.35, 0.65, 0.95),
    `5` = c(0.05, 0.275, 0.5, 0.725, 0.95),
    `6` = c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95),
    `7` = c(0.025, 0.1833, 0.3417, 0.5, 0.6583, 0.8167, 0.975)
  )
  as.numeric(stats::quantile(x, pr, names = FALSE))
}

# The background of Childs et al. (2022, equation 2): for each row, the
# median PM2.5 of its place and calendar month over the years y - 1 to y +
# 1 on days without smoke; NA (with a warning) when there are none.
.axchilds <- function(d, unit, day, pv, pl, fn) {
  if (!inherits(d[[day]], "Date")) {
    .ckstop(fn, "`cf = \"childs\"` needs `day` as dates")
  }
  yr <- as.integer(format(d[[day]], "%Y"))
  um <- paste(as.character(d[[unit]]), format(d[[day]], "%m"), sep = "\r")
  bg <- rep(NA_real_, length(pv))
  for (g in split(seq_along(pv), um)) {
    for (y in unique(yr[g])) {
      w <- g[abs(yr[g] - y) <= 1L & pl[g] == 0]
      if (length(w)) bg[g[yr[g] == y]] <- stats::median(pv[w])
    }
  }
  if (anyNA(bg)) {
    warning(fn, ": ", sum(is.na(bg)), " row(s) have no smoke-free days ",
            "in their place, month and three years; their background is NA.",
            call. = FALSE)
  }
  bg
}

# The season of each row for smkiv()'s backgrounds: the column seas, or
# the meteorological season of a date day, or one season.
.axsmsea <- function(d, day, seas, pm, f) {
  if (!is.null(seas)) {
    as.character(d[[seas]])
  } else {
    if (!inherits(d[[day]], "Date") && !is.null(pm)) {
      .ckstop(
        f,
        "give `seas`, or `day` as dates for the meteorological seasons"
      )
    }
    if (inherits(d[[day]], "Date")) {
      m <- as.integer(format(d[[day]], "%m"))
      c("DJF", "MAM", "JJA", "SON")[(m %/% 3) %% 4 + 1L]
    } else {
      rep("all", nrow(d))
    }
  }
}

# The PM2.5 columns of smkiv(): the background of each place and season
# from its days without a plume, the indicator adjusted by its own or its
# neighbour's PM2.5, and the smoke PM2.5 above the background.
.axsmpm <- function(d, out, pl, key, pm, adj, nbr, z, cf, unit, day, f) {
  pv <- .axmat(d, pm, f)[, 1L]
  cl <- pl == 0
  if (any(tapply(cl, key, sum) < 2L)) {
    .ckstop(
      f,
      paste0(
        "a place and season has fewer than two days without a plume for ",
        "its background"
      )
    )
  }
  mu <- tapply(pv[cl], key[cl], mean)
  sdv <- tapply(pv[cl], key[cl], stats::sd)
  md <- tapply(pv[cl], key[cl], stats::median)
  bmu <- as.numeric(mu[key])
  bsd <- as.numeric(sdv[key])
  hi <- pl == 1 & pv > bmu + z * bsd
  if (adj == "own") out$smoke <- as.numeric(hi)
  if (adj == "nbr") {
    if (is.null(nbr)) {
      .ckstop(f, "`adj = \"nbr\"` needs `nbr`, the neighbour of each row")
    }
    .axcols(d, nbr, "nbr", f, n = 1L)
    j <- match(
      paste(as.character(d[[nbr]]), as.character(d[[day]]), sep = "\r"),
      paste(as.character(d[[unit]]), as.character(d[[day]]), sep = "\r")
    )
    out$smoke <- as.numeric(hi[j])
  }
  if (cf == "childs") {
    out$bg <- .axchilds(d, unit, day, pv, pl, f)
    out$spm <- ifelse(pl == 1, pmax(pv - out$bg, 0), 0)
  } else {
    out$bg <- if (cf == "mean") bmu else as.numeric(md[key])
    out$spm <- ifelse(!is.na(out$smoke) & out$smoke == 1,
                      pmax(pv - out$bg, 0), 0)
  }
  out
}

# smkiv()'s columns aggregated to the groups of by: means of the
# PM2.5 columns, sums of the others.
.axsmby <- function(out, d, by) {
  g <- .axcode(d[[by]])
  a <- data.frame(by = unique(as.character(d[[by]])))
  names(a) <- by
  for (k in names(out)) {
    v <- out[[k]]
    a[[k]] <- if (k %in% c("bg", "spm")) {
      as.numeric(tapply(v, g$c, mean))
    } else {
      as.numeric(tapply(v, g$c, sum))
    }
  }
  a
}

# hetfx(how = "split"): est on the rows of each level of the moderator:
# list(tab, b the estimates, v their covariance, diagonal).
.axhtspl <- function(d, est, a, var, mv, lv, f) {
  rows <- lapply(lv, function(l) {
    r <- do.call(est, c(list(data = d[mv == l, , drop = FALSE]), a))
    cf <- .axrow(r, var, f)
    data.frame(
      level = l,
      est = cf$est,
      se = cf$se,
      lo = cf$lo,
      hi = cf$hi,
      p = cf$p,
      n = r$n
    )
  })
  tab <- do.call(rbind, rows)
  b <- tab$est
  v <- diag(tab$se^2, length(b))
  list(tab = tab, b = b, v = v)
}

# hetfx(how = "inter"): one fit of est with the exposure (and the
# instruments) times the indicator of each level of the moderator:
# list(tab, b the estimates, v their covariance).
.axhtint <- function(d, est, a, var, mv, lv) {
  nm <- paste0(var, ":", lv)
  for (j in seq_along(lv)) d[[nm[j]]] <- d[[var]] * (mv == lv[j])
  zi <- NULL
  if (!is.null(a[["inst"]])) {
    zi <- character()
    for (z in a[["inst"]]) {
      for (j in seq_along(lv)) {
        zn <- paste0(z, ":", lv[j])
        d[[zn]] <- d[[z]] * (mv == lv[j])
        zi <- c(zi, zn)
      }
    }
  }
  a <- .axrepl(a, var, nm, zi)
  r <- do.call(est, c(list(data = d), a))
  vc <- .axvcov(r)
  cf <- r$coef[match(nm, r$coef$term), , drop = FALSE]
  tab <- data.frame(
    level = lv,
    est = cf$est,
    se = cf$se,
    lo = cf$lo,
    hi = cf$hi,
    p = cf$p,
    n = as.integer(table(factor(mv, levels = lv)))
  )
  b <- cf$est
  v <- vc[nm, nm, drop = FALSE]
  list(tab = tab, b = b, v = v)
}

# The exposures of ccross(): the columns of lagv lag days before, and the
# columns of ma averaged over k days, a place's day left out where one of
# them is missing.
.axccma <- function(d, unit, day, ma, k, lagv, lag) {
  if (length(lagv) && lag > 0) {
    # the value lag days before, of the same place (a day without one is
    # dropped)
    ix <- .Call(C_clab_r_dshift, .axcode(d[[unit]])$c,
                as.numeric(d[[day]]), -lag)
    for (v in lagv) d[[v]] <- d[[v]][ix]
    d <- d[stats::complete.cases(d[lagv]), , drop = FALSE]
  }
  if (length(ma) && k > 1) {
    uc <- .axcode(d[[unit]])$c
    tv <- as.numeric(d[[day]])
    s <- lapply(ma, function(v) d[[v]])
    for (j in seq_len(k - 1)) {
      ix <- .Call(C_clab_r_dshift, uc, tv, -j)
      for (m in seq_along(ma)) s[[m]] <- s[[m]] + d[[ma[m]]][ix]
    }
    for (m in seq_along(ma)) d[[ma[m]]] <- s[[m]] / k
    d <- d[stats::complete.cases(d[ma]), , drop = FALSE]
  }
  d
}

# The natural cubic splines of ccross() of the covariates in ns: with d
# degrees of freedom, knots at the range and the quantiles 1/d, ...,
# (d-1)/d (those of splines::ns(x, df = d)), in the restricted cubic
# spline basis, which spans the same functions: list(d with the basis
# columns, a with them among the regressors x).
.axccns <- function(d, a, ns, f) {
  for (v in names(ns)) {
    x <- d[[v]]
    if (!is.numeric(x)) .ckstop(f, "column `%s` must be numeric", v)
    dg <- ns[[v]]
    kn <- c(min(x), stats::quantile(x, seq_len(dg - 1L) / dg, names = FALSE),
            max(x))
    b <- .axrcs(x, kn)
    nm <- paste0(v, ":ns", seq_len(ncol(b)))
    for (j in seq_len(ncol(b))) d[[nm[j]]] <- b[, j]
    a[["x"]] <- c(a[["x"]], nm)
  }
  list(d = d, a = a)
}

# The checks of jntx()'s columns and numbers, in their order:
# list(k1, k2, p1, p2, ngrid).
.axjxck <- function(d, y, x1, x2, est, k1, k2, p1, p2, ngrid, f) {
  .axcols(d, c(y, x1, x2), "y, x1, x2", f, n = 3L)
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  k1 <- .cknum(k1, "k1", f, min = 3, max = 7, whole = TRUE)
  k2 <- .cknum(k2, "k2", f, min = 3, max = 7, whole = TRUE)
  p1 <- .cknum(p1, "p1", f, scalar = FALSE, min = 0, max = 1)
  p2 <- .cknum(p2, "p2", f, scalar = FALSE, min = 0, max = 1)
  if (length(p1) != 2L) .ckstop(f, "`p1` must be two quantiles")
  if (length(p2) < 2L) .ckstop(f, "`p2` must be two or more quantiles")
  ngrid <- .cknum(ngrid, "ngrid", f, min = 2, whole = TRUE)
  if (ngrid > 2^13) {
    .ckstop(f, "`ngrid` times `ngrid` must be at most 67108864")
  }
  list(k1 = k1, k2 = k2, p1 = p1, p2 = p2, ngrid = ngrid)
}

# The basis of jntx() at u and v: the two exposures and their product
# (linear), or their restricted cubic splines at the knots t1 and t2 and
# the products of the two bases (spline).
.axjxbas <- function(u, v, form, t1, t2) {
  if (form == "linear") return(cbind(u, v, u * v))
  b1 <- .axrcs(u, t1)
  b2 <- .axrcs(v, t2)
  cbind(
    b1,
    b2,
    do.call(cbind, lapply(seq_len(ncol(b2)), function(j) b1 * b2[, j]))
  )
}

# What jntx() reports from rel(), the change of the linear predictor
# between two points: the effect of x1 from its p1 quantiles at the p2
# quantiles of x2 (cond), and over a grid of both against ref (surf).
.axjxout <- function(a, p1, p2, ref, ngrid, rel) {
  q1 <- as.numeric(stats::quantile(a[, 1L], p1, names = FALSE))
  q2 <- as.numeric(stats::quantile(a[, 2L], p2, names = FALSE))
  cc <- rel(rep(q1[2L], length(q2)), q2, rep(q1[1L], length(q2)), q2)
  z <- stats::qnorm(0.975)
  cond <- data.frame(
    x2 = q2,
    p2 = p2,
    est = cc$e,
    se = cc$s,
    lo = cc$e - z * cc$s,
    hi = cc$e + z * cc$s,
    rr = exp(cc$e)
  )
  g <- expand.grid(x1 = seq(min(a[, 1L]), max(a[, 1L]), length.out = ngrid),
                   x2 = seq(min(a[, 2L]), max(a[, 2L]), length.out = ngrid))
  sf <- rel(g$x1, g$x2, rep(ref[1L], nrow(g)), rep(ref[2L], nrow(g)))
  g$est <- sf$e
  g$se <- sf$s
  g$rr <- exp(sf$e)
  list(cond = cond, surf = g)
}

# The bins of dosebin(): nq quantile bins of xv or the breaks brk, two at
# least, every value in one: list(bn the bin of each row, lv, nbin).
.axdbbin <- function(xv, nq, brk, var, f) {
  if (is.null(brk)) {
    nq <- .cknum(nq, "nq", f, min = 2, whole = TRUE)
    if (nq > 2^26) .ckstop(f, "`nq` must be at most 67108864")
    qs <- seq(0, 1, length.out = nq + 1L)
    brk <- unique(as.numeric(stats::quantile(xv, qs, names = FALSE)))
    # cut() reads one number as a count of intervals, not as a break: an
    # exposure with one value (or whose quantiles all fall on one) has no
    # bins to compare
    if (length(brk) < 3L) {
      .ckstop(f, "`%s` gives fewer than two bins", var)
    }
  } else if (!is.numeric(brk) || length(brk) < 3L || anyNA(brk) ||
               is.unsorted(brk, strictly = TRUE)) {
    .ckstop(f, "`brk` must be three or more increasing numbers")
  }
  bn <- cut(xv, brk, include.lowest = TRUE)
  if (anyNA(bn)) .ckstop(f, "some values of `%s` fall outside `brk`", var)
  lv <- levels(bn)
  nbin <- length(lv)
  if (nbin < 2L) .ckstop(f, "`%s` gives fewer than two bins", var)
  list(bn = bn, lv = lv, nbin = nbin)
}

# The Wald tests of dosebin(): every bin equal to the reference (all),
# and, with two bins or more, the bins on the line through the reference
# and the first bin at their medians (linear).
.axdbtst <- function(b, v, mid, use, ref, k) {
  tst <- .axwald(b, v, diag(k))
  tst$test <- "all"
  if (k >= 2L) {
    dd <- mid[use] - mid[ref]
    rl <- t(vapply(2:k, function(q) {
      w <- numeric(k)
      w[q] <- 1 / dd[q]
      w[1L] <- -1 / dd[1L]
      w
    }, numeric(k)))
    tl <- .axwald(b, v, rl)
    tl$test <- "linear"
    tst <- rbind(tst, tl)
  }
  tst
}

# The checks of smkiv()'s columns and numbers, in their order: list(pl
# the plume indicator, z).
.axsmck <- function(d, plume, unit, day, pm, seas, by, dens, adj, z, f) {
  .axcols(d, c(plume, unit, day), "plume, unit, day", f, n = 3L)
  .axcols(d, pm, "pm", f, n = 1L, null = TRUE)
  .axcols(d, seas, "seas", f, n = 1L, null = TRUE)
  .axcols(d, by, "by", f, n = 1L, null = TRUE)
  if (!is.null(dens)) {
    if (!is.character(dens) || length(dens) != 1L || !dens %in% names(d)) {
      .ckstop(f, "`dens` must name a column of `data`")
    }
  }
  pl <- .axmat(d, plume, f)[, 1L]
  if (!all(pl %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", plume)
  z <- .cknum(z, "z", f, min = 0)
  if (adj != "raw" && is.null(pm)) .ckstop(f, "`adj = \"%s\"` needs `pm`", adj)
  list(pl = pl, z = z)
}

# smkiv()'s indicators of the levels of the density column dens, added
# to out (out unchanged without dens).
.axsmdn <- function(out, d, dens) {
  if (!is.null(dens)) {
    dv <- as.character(d[[dens]])
    lv <- sort(unique(dv[!is.na(dv) & dv != ""]))
    for (l in lv) out[[paste0("dens:", l)]] <- as.numeric(!is.na(dv) & dv == l)
  }
  out
}

# The checks of ccross()'s columns, estimator and numbers, in their
# order: list(k, lag).
.axccck <- function(d, y, unit, day, ma, est, k, lagv, lag, ns, f) {
  .axcols(d, c(y, unit, day), "y, unit, day", f, n = 3L)
  .axcols(d, ma, "ma", f, null = TRUE)
  if (!inherits(d[[day]], "Date")) {
    .ckstop(f, "column `%s` must hold dates", day)
  }
  if (!.axest(est)) {
    .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  }
  k <- .cknum(k, "k", f, min = 1, whole = TRUE)
  .axcols(d, lagv, "lagv", f, null = TRUE)
  lag <- .cknum(lag, "lag", f, min = 0, whole = TRUE)
  if (!is.null(ns)) {
    if (!is.numeric(ns) || is.null(names(ns)) || anyNA(ns) ||
          any(ns < 1 | ns != round(ns))) {
      .ckstop(f, "`ns` must be whole degrees of freedom named by column")
    }
    .axcols(d, names(ns), "names(ns)", f)
  }
  list(k = k, lag = lag)
}
