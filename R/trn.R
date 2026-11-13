#' Trend in a series: rate ratio, Mann-Kendall and Theil-Sen
#'
#' Measures the trend in a short series of counts (admissions, incidents
#' or placements per year) three ways. The Poisson log-linear model
#' \eqn{\log \mu_t = b_0 + b_1 (x_t - \bar x) + \log e_t} gives the rate
#' ratio per unit of `x`, `exp(b1)` (0.85: a fall of 15% a period), with
#' the exposure `expo` as offset so that the trend is in the rate; when the
#' Pearson dispersion (over `n - 2`) is above 1.5 the interval and p-value
#' are the quasi-Poisson ones (standard error times the square root of the
#' dispersion, the t distribution on `n - 2` degrees of freedom). The fit
#' takes the steps of [stats::glm()]. Mann-Kendall's rank test and the
#' Theil-Sen slope (the median of the pairwise slopes, with Sen's
#' distribution-free interval) make no distributional assumption; the
#' p-value is exact over every ordering of the observed values up to 8
#' periods (ties kept) and normal with the tie and continuity corrections
#' beyond. Negating the series negates S, so rises and falls are treated
#' alike (`proof("trndrr")`). `yoy` holds the change on the period before,
#' and with `seas` the
#' totals by season are tested for equality (chi-square).
#'
#' When every count falls in the first or the last period the likelihood
#' rises without bound: the rate ratio is 0 or `Inf`, with `NA` for the
#' rest of its row, and a warning. Sen's interval needs enough periods:
#' `clamp` is `TRUE` when an end had to be the extreme slope.
#'
#' @param y The counts, in period order: numbers from 0, not all 0, at
#'   least 3.
#' @param x `NULL` for the periods 1, 2, ..., or the periods: strictly
#'   increasing numbers (years, say).
#' @param expo `NULL` or the exposure of each period (population,
#'   person-days): positive numbers.
#' @param seas `NULL` or the season of each count (month, day of week):
#'   labels, two or more.
#' @param alt The alternative of the Mann-Kendall test: `"both"`, `"up"`
#'   or `"down"`.
#' @inheritParams kmrec
#' @return An object of class `trndrr`, a list: `rr`, the rate ratio row
#'   (`rr`, `lo`, `hi`, `p`, the log slope `b` and its `se`, the dispersion
#'   `disp`, `over` when the quasi-Poisson interval was used); `mk`
#'   (`s`, its variance `var`, `tau`, `p`, `exact`); `sen` (`slope`, `int`,
#'   `lo`, `hi`, `clamp`); `yoy` (`x`, `y`, the change `chg` and per cent
#'   change `pct` on the period before); `seas`, `NULL` or a list of `tab`
#'   (`season`, `n`, `expd`), `stat`, `df` and `p`; `fitted`, the Poisson
#'   fitted counts; `n`, `alt` and `level`.
#' @references
#'   Sen, P. K. (1968). Estimates of the regression coefficient based on
#'   Kendall's tau. *Journal of the American Statistical Association*,
#'   63(324), 1379-1389. \doi{10.1080/01621459.1968.10480934}
#'
#'   Hamed, K. H. and Rao, A. R. (1998). A modified Mann-Kendall trend test
#'   for autocorrelated data. *Journal of Hydrology*, 204(1-4), 182-196.
#'   \doi{10.1016/S0022-1694(97)00125-X}
#' @seealso [stpcp()] for a step change, [xcser()] to correlate two series.
#' @examples
#' y <- c(402, 377, 390, 368, 331, 310, 305)
#' trndrr(y, 2017:2023)
#' # a rate: counts over a growing population
#' trndrr(c(20, 25, 30, 33), expo = c(1000, 1500, 2500, 2600))$rr
#' # equal totals by quarter?
#' trndrr(c(5, 9, 7, 4, 6, 10, 8, 5), seas = rep(1:4, 2))$seas
#' @export
trndrr <- function(
  y,
  x = NULL,
  expo = NULL,
  seas = NULL,
  alt = "both",
  level = 0.95
) {
  .Call(C_clab_r_trndrr, y, x, expo, seas, alt, level)
}

#' @export
print.trndrr <- function(x, ...) {
  r <- x[["rr"]]
  cat("trndrr: ", x[["n"]], " periods\nrate ratio per period ",
      format(r[["rr"]]), " (", format(100 * x[["level"]]), "% ",
      format(r[["lo"]]), " to ", format(r[["hi"]]), "), p ", format(r[["p"]]),
      if (isTRUE(r[["over"]])) "; overdispersed: quasi-Poisson interval",
      "\nMann-Kendall S ", format(x[["mk"]][["s"]]), ", p ",
      format(x[["mk"]][["p"]]), "; Theil-Sen slope ",
      format(x[["sen"]][["slope"]]), "\n", sep = "")
  if (!is.null(x[["seas"]])) {
    cat("seasons: chi-square ", format(x[["seas"]][["stat"]]), " on ",
        x[["seas"]][["df"]], " df, p ", format(x[["seas"]][["p"]]), "\n",
        sep = "")
  }
  invisible(x)
}

#' Change point: a step in a series and Pettitt's test
#'
#' Asks whether a series moved to a new level, and where. The step scan
#' takes every split after period `k` with at least `minseg` periods on
#' each side, and the largest between-segment sum of squares \eqn{T_k = k
#' (n - k) / n (\bar y_a - \bar y_b)^2}; its p-value comes from the same
#' maximum over reorderings of the series (every ordering when `n!` is at
#' most `nperm`, giving an exact p-value; otherwise `nperm` random ones,
#' `(1 + #{T* >= T}) / (1 + nperm)`), so it is a test of a step anywhere,
#' not of the best split alone. Pettitt's rank test takes the split with
#' the largest \eqn{|U_t|}, \eqn{U_t = \sum_{i \le t < j} sign(y_j -
#' y_i)}, with its approximate p-value \eqn{2 \exp(-6 K^2 / (n^3 + n^2))}.
#' Both assume exchangeable periods under no change: a steady trend can
#' read as a step.
#'
#' @param y The series, in period order: finite numbers.
#' @param x `NULL`, or the periods (any vector, one per value): `at` then
#'   names the last period before each change.
#' @param minseg The fewest periods on either side of a step: one whole
#'   number from 1.
#' @param nperm The reorderings of the step scan's null: one whole number
#'   from 1.
#' @inheritParams hwksim
#' @return An object of class `stpcp`, a list: `tab`, one row for each
#'   test (`test`, `k`, the last period before the change, the means
#'   `before` and `after`, the statistic `stat` (T or K) and `p`, with `at`
#'   from `x`); `nperm`, the orderings used; `exact`; and `n`.
#' @references
#'   Pettitt, A. N. (1979). A non-parametric approach to the change-point
#'   problem. *Journal of the Royal Statistical Society, Series C*, 28(2),
#'   126-135. \doi{10.2307/2346729}
#' @examples
#' stpcp(c(100, 104, 98, 60, 63, 58), x = 2018:2023)
#' stpcp(c(100, 104, 98, 101, 97, 60, 63, 58, 61, 59), nperm = 999, seed = 1)
#' @export
stpcp <- function(y, x = NULL, minseg = 2, nperm = 9999, seed = NULL) {
  r <- .Call(C_clab_r_stpcp, y, minseg, nperm, seed)
  if (!is.null(x)) {
    if (length(x) != r[["n"]]) {
      stop("stpcp: `x` must have one period per value, ", r[["n"]], ".",
           call. = FALSE)
    }
    r[["tab"]][["at"]] <- x[r[["tab"]][["k"]]]
  }
  r
}

#' @export
print.stpcp <- function(x, ...) {
  cat("stpcp: ", x[["n"]], " periods; step p-value ",
      if (x[["exact"]]) "exact over " else "from ", format(x[["nperm"]]),
      " orderings\n", sep = "")
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}

#' CUSUM surveillance of counts
#'
#' Watches a series of counts for a rise above what is expected: the
#' cumulative sum \eqn{s_t = \max(0, s_{t-1} + z_t - k)} of the
#' standardised excess \eqn{z_t} raises an alarm where \eqn{s_t \ge h}.
#' The expected counts are given (`expd`), or are the mean of the counts
#' before `start`, the first period watched. `z` is `(y - m) / sqrt(m)`
#' (`"std"`), Rossi's transform `(y - 3 m + 2 sqrt(y m)) / (2 sqrt(m))`
#' (`"rossi"`), Anscombe's `1.5 (y^(2/3) - m^(2/3)) / m^(1/6)`
#' (`"anscombe"`) or the count itself (`"none"`, for `k` and `h` on the
#' count scale). With `reset` the sum starts again from 0 after an alarm.
#'
#' @param y The counts: numbers from 0.
#' @param expd `NULL`, or the expected count of each period: positive
#'   numbers.
#' @param start Without `expd`, the first period watched: one whole number
#'   from 2.
#' @param k The reference value subtracted each period: one number.
#' @param h The decision limit: one positive number.
#' @param trans The standardisation: `"std"`, `"rossi"`, `"anscombe"` or
#'   `"none"`.
#' @param reset `TRUE` to start again after each alarm.
#' @return A data frame, one row per period watched: `i`, its position;
#'   `y`; `expd`; `z`; the sum `s`; and `alarm`.
#' @references
#'   Page, E. S. (1954). Continuous inspection schemes. *Biometrika*,
#'   41(1-2), 100-115. \doi{10.1093/biomet/41.1-2.100}
#'
#'   Rossi, G., Lampugnani, L. and Marchi, M. (1999). An approximate CUSUM
#'   procedure for surveillance of health events. *Statistics in Medicine*,
#'   18(16), 2111-2122.
#' @examples
#' y <- c(10, 12, 9, 11, 10, 13, 19, 22, 25, 12)
#' cusrv(y, start = 6)
#' cusrv(y, expd = rep(11, 10), trans = "rossi", reset = TRUE)
#' @export
cusrv <- function(
  y,
  expd = NULL,
  start = NULL,
  k = 1.04,
  h = 2.26,
  trans = "std",
  reset = FALSE
) {
  .Call(C_clab_r_cusrv, y, expd, start, k, h, trans, reset)
}

#' Correlation of two series, at lags
#'
#' Correlates two series period by period: Pearson's r with Fisher's
#' interval and the t-test p-value of [stats::cor.test()], and the
#' cross-correlations at lags \eqn{-K..K} as [stats::ccf()] estimates them
#' (lag `k` is the correlation of `x[t + k]` with `y[t]`; the default `K`
#' is `ccf()`'s, `floor(10 log10(n / 2))`). Two named series are matched by
#' name (the periods both have, in the order of `x`), so counts from two
#' sources for the same years line up whatever years each covers. With
#' `dif` the changes from one period to the next are correlated instead of
#' the levels: two series that both trend correlate in level whether or
#' not they move together. Every correlation lies in \eqn{[-1, 1]}
#' (Cauchy-Schwarz, `proof("xcser")`).
#'
#' @param x,y The series: finite numbers, of one length, or both named by
#'   period.
#' @param lagmax `NULL` or the largest lag: one whole number from 0, below
#'   the number of periods.
#' @param dif `TRUE` to correlate the period-to-period changes.
#' @inheritParams kmrec
#' @return An object of class `xcser`, a list: `cor` (`r`, `lo`, `hi`, `p`,
#'   `n`); `ccf` (`lag`, `r`); `common`, the periods matched by name
#'   (`NULL` without names); and `level`.
#' @examples
#' a <- c("2019" = 120, "2020" = 98, "2021" = 105, "2022" = 130, "2023" = 141)
#' b <- c("2018" = 50, "2019" = 61, "2020" = 44, "2021" = 52, "2022" = 66,
#'        "2023" = 70)
#' xcser(a, b)
#' xcser(a, b, dif = TRUE)$cor
#' @export
xcser <- function(x, y, lagmax = NULL, dif = FALSE, level = 0.95) {
  .Call(C_clab_r_xcser, x, y, lagmax, dif, level)
}

#' @export
print.xcser <- function(x, ...) {
  r <- x[["cor"]]
  cat("xcser: ", r[["n"]], " periods; r ", format(r[["r"]]), " (",
      format(100 * x[["level"]]), "% ", format(r[["lo"]]), " to ",
      format(r[["hi"]]), "), p ", format(r[["p"]]), "\n", sep = "")
  print(x[["ccf"]], row.names = FALSE)
  invisible(x)
}
