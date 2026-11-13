#' An event study: the effect period by period around adoption
#'
#' Leads and lags of a policy's adoption from a two-way fixed-effects
#' regression on indicators of the time since adoption, with the joint test
#' that the leads are zero (no effect before the policy, the testable part
#' of parallel trends), and the raw mean outcome of each adoption cohort
#' around its adoption.
#'
#' With `G_i` the period unit `i` is first treated and `k = t - G_i`, the
#' model is `y_it = sum_k g_k 1(t - G_i = k) + adj_it'th + a_i + m_t + e_it`
#' for `k` from `-leads` to `lags` but `ref`, the period left out (the one
#' before adoption by default). Times beyond the window are binned into its
#' ends (`bin = TRUE`) or their rows dropped. Units never treated have every
#' indicator 0. The variance, the clusters and the degrees of freedom are
#' those of [didfe()]. The pre-trend test is the Wald statistic `W =
#' g_pre' V_pre^-1 g_pre` over the leads, reported as `F = W / q` on `q` and
#' `G - 1` degrees of freedom.
#'
#' With staggered adoption and effects that change over time, the lead and
#' lag coefficients can mix in effects from other periods (Sun and Abraham
#' 2021); [stgdid()] gives event-time effects free of that.
#'
#' `coh` is the mean outcome of each cohort (the units first treated in one
#' period) in each period, as a trailing moving average over `ma` periods,
#' by time since adoption within the window, and with `norm = TRUE` scaled
#' to 1 in the period before adoption: a check, before any model, that the
#' cohorts moved together before the policy.
#'
#' With `mbar`, the effects are bounded when the trend after adoption may
#' depart from parallel by up to `mbar` times the largest departure seen
#' before it (the relative magnitudes of Rambachan and Roth 2023): with
#' `D` the largest change between consecutive pre-period coefficients (the
#' reference at 0), the bias at event time `e >= 0` is at most `(e + 1)
#' mbar D`, so the identified set of the target (the mean effect after
#' adoption, or the event time `target`) is its estimate `+-` the mean of
#' `(e + 1)` times `mbar D` (`lob`, `upb`). The interval (`lo`, `hi`) draws
#' the coefficients `draws` times from their normal distribution with the
#' cluster-robust covariance, computes the bounds for each, and takes the
#' `(1 - level) / 2` quantile of the lower and the `(1 + level) / 2`
#' quantile of the upper: by Bonferroni it covers the identified set with
#' probability at least `level`. It is conservative; the conditional and
#' hybrid intervals of Rambachan and Roth are shorter. `brk` is the
#' smallest `mbar` (by bisection, up to 10; `Inf` beyond) whose interval
#' holds 0.
#'
#' @inheritParams didfe
#' @param time The period of each row: whole numbers (a year, a month
#'   counted from a start); each unit at most once per period.
#' @param first The period each row's unit is first treated: numbers,
#'   constant within units; `NA`, 0 or `Inf` for a unit never treated.
#'   Treatment is taken to stay on.
#' @param leads,lags The periods before and after adoption estimated: whole
#'   numbers, `leads` at least 1 and `lags` at least 0.
#' @param ref The event time left out: a whole number from `-leads` to
#'   `lags`.
#' @param bin `TRUE` to bin event times beyond the window into its ends,
#'   `FALSE` to drop their rows.
#' @param ma The width of the moving average of the cohort means: a whole
#'   number, at least 1 (1 for none).
#' @param norm `TRUE` to scale each cohort's means to 1 in the period
#'   before adoption.
#' @param mbar `NULL`, or the bounds on the departure from parallel trends
#'   relative to the largest before adoption: numbers from 0.
#' @param target `NULL` for the mean effect after adoption, or one event
#'   time from 0.
#' @param draws The normal draws of the sensitivity interval: a whole
#'   number from 100.
#' @param seed The seed of the draws: `NULL` to draw it from R's generator,
#'   or a whole number from 0.
#' @return An object of class `evtst`, a list: `est`, one row per event
#'   time (`evtm`, `est`, `se`, `lo`, `hi`, `stat`, `p`), the reference a
#'   row with `est` 0; `pre`, the pre-trend test (`f`, `df1`, `df2`, `p`;
#'   `NA` without leads to test); `vcov`, the variance of the estimated
#'   event times; `coh`, the cohort means (`cohort`, `evtm`, `value`); `n`,
#'   the rows used; `ncl`, the clusters; and `level`; with `mbar`, `sens`
#'   (`mbar`, `lob`, `upb`, `lo`, `hi`) and `brk`, else `NULL`.
#' @references
#'   Sun, L. and Abraham, S. (2021). Estimating dynamic treatment effects in
#'   event studies with heterogeneous treatment effects. *Journal of
#'   Econometrics*, 225(2), 175-199. \doi{10.1016/j.jeconom.2020.09.006}
#'
#'   Freyaldenhoven, S., Hansen, C. and Shapiro, J. M. (2019). Pre-event
#'   trends in the panel event-study design. *American Economic Review*,
#'   109(9), 3307-3338. \doi{10.1257/aer.20180609}
#' @examples
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
#' i <- match(d$town, LETTERS)
#' d$first <- c(2019, 2019, 2019, 2019, 2021, 2021, 2021, 2021,
#'              NA, NA, NA, NA)[i]
#' on <- !is.na(d$first) & d$year >= d$first
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * on + sin(seq_len(120))
#' e <- evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3)
#' e
#' e$pre
#' # how far trends may bend before the effect is lost
#' evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3,
#'       mbar = c(0, 1, 2), seed = 1)$sens
#' @export
evtst <- function(y, unit, time, first, leads = 4, lags = 4, ref = -1,
                  bin = TRUE, adj = NULL, clus = NULL, w = NULL, ma = 1,
                  norm = TRUE, mbar = NULL, target = NULL, draws = 5000,
                  seed = NULL, level = 0.95) {
  .Call(C_clab_r_evtst, y, unit, time, first, leads, lags, ref, bin, adj,
        clus, w, ma, norm, mbar, target, draws, seed, level)
}

#' @export
print.evtst <- function(x, ...) {
  cat("evtst: ", x[["n"]], " rows, ", x[["ncl"]], " clusters\n", sep = "")
  print(x[["est"]], row.names = FALSE)
  p <- x[["pre"]]
  if (!is.na(p[["f"]])) {
    cat("leads jointly zero: F = ", format(p[["f"]]), " on ", p[["df1"]],
        " and ", p[["df2"]], " df, p = ", format(p[["p"]]), "\n", sep = "")
  }
  if (!is.null(x[["sens"]])) {
    cat("bounds under departures from parallel trends:\n")
    print(x[["sens"]], row.names = FALSE)
    cat("breakdown mbar ", format(x[["brk"]]), "\n", sep = "")
  }
  invisible(x)
}
