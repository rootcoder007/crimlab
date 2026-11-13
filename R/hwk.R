#' Hawkes process: self-exciting event times
#'
#' Fits a self-exciting (Hawkes) point process to event times, such as the
#' times of shootings, burglaries or calls for service, by maximum
#' likelihood: a background rate and, on top of it, a rise in the rate
#' after each event that dies away with the time since. The fitted
#' branching ratio `eta` is the expected number of events each event
#' triggers; it says how much of the record is contagion rather than
#' background.
#'
#' The intensity is \deqn{\lambda(t) = \nu(t) + \eta \sum_{t_j < t} g(t -
#' t_j)} over the window `[0, hor]`, with a normalised kernel `g`:
#' exponential (`"exp"`, rate `beta`: a Markovian process), Weibull
#' (`"weibull"`, `shape` and `scale`) or Lomax (`"lomax"`, the power-law
#' kernel of seismology's modified Omori law, `shape` and `scale`, density
#' `shape / scale (1 + t / scale)^-(shape + 1)`). The baseline is constant
#' (`"const"`, \eqn{\nu = e^{a_0}}) or a log-linear trend with one cycle of
#' length `per` (`"cycle"`, \eqn{\nu(t) = \exp(a_0 + a_1 t / hor + a_2
#' \sin(2 \pi t / per) + a_3 \cos(2 \pi t / per))}), whose integral is taken
#' by Simpson's rule on steps of at most `per / 16`. The log-likelihood is
#' \eqn{\sum_i \log \lambda(t_i) - \int_0^{hor} \lambda}; it is maximised
#' with its analytic gradient by BFGS over the baseline terms, the logit of
#' `eta` (so `0 < eta < 1`) and the logarithms of the kernel parameters.
#' The exponential kernel is fitted from two starts and the better optimum
#' kept; the others from their own start and from the exponential fit (the
#' shape-1 member). Standard errors come from the inverse of the numerical
#' Hessian of the analytic gradient.
#'
#' Several kernels or baselines fit every combination: `cmp` compares them
#' (the exponential kernel against the others is the comparison of a
#' Markovian with a non-Markovian process) and the one with the least AIC
#' is returned. `branch` holds the arithmetic of the branching process:
#' a cluster started by one background event has `1 / (1 - eta)` events on
#' average, a constant baseline \eqn{\nu} gives the stationary rate
#' \eqn{\nu / (1 - eta)}, and the share of events that are triggered is
#' `eta` (`proof("hwkfit")`).
#'
#' Events at one instant do not excite each other: the sum runs over
#' strictly earlier events. Times recorded to a resolution (dates, whole
#' hours) carry ties the process never makes; `res` spreads each time
#' uniformly over `[t, t + res)` first (cut at `hor`), from `seed`.
#' The exponential kernel costs O(n) per likelihood, the others O(n^2).
#'
#' @param times Event times: numbers in one unit (days, say), sorted, in
#'   `[0, hor]`; at least 10 of them.
#' @param hor The end of the observation window: one positive number.
#' @param kernel The triggering kernel: one or more of `"exp"`,
#'   `"weibull"` and `"lomax"`.
#' @param base The baseline: one or both of `"const"` and `"cycle"`.
#' @param per The length of the cycle of the `"cycle"` baseline, in the
#'   unit of `times` (365.25 for a yearly cycle in days).
#' @param start `NULL`, or the starting parameters of a single kernel and
#'   baseline in the order of `par` (`eta` between 0 and 1).
#' @param res `NULL`, or the resolution the times were recorded to: one
#'   positive number.
#' @param seed `NULL` to draw the seed of `res`'s spreading from R's
#'   generator, or one whole number from 0.
#' @param maxit The most BFGS steps per start: one whole number from 1.
#' @return An object of class `hwkfit`, a list: `par`, the estimates
#'   (`a0` or `a0` to `a3` on the log scale, `eta`, then `beta` or `shape`
#'   and `scale`); `se` and `cov`, their standard errors and covariance;
#'   `loglik`, `aic`, `bic`; `branch` (`eta`, the mean cluster `size`, the
#'   stationary `rate` of a constant baseline and the triggered `share`);
#'   `conv` (`TRUE` when the gradient vanished) and `iter`; `kernel`,
#'   `base`, `per`, `n`, `hor`, `jit` (whether `res` spread the times) and
#'   `times` (the times fitted); and `cmp`, one row per kernel and baseline
#'   fitted (`kernel`, `base`, `loglik`, `k` parameters, `aic`, `bic`,
#'   `eta`, `conv`). A branching ratio above 0.99 is warned about.
#' @references
#'   Hawkes, A. G. (1971). Spectra of some self-exciting and mutually
#'   exciting point processes. *Biometrika*, 58(1), 83-90.
#'   \doi{10.1093/biomet/58.1.83}
#'
#'   Ozaki, T. (1979). Maximum likelihood estimation of Hawkes'
#'   self-exciting point processes. *Annals of the Institute of Statistical
#'   Mathematics*, 31, 145-155. \doi{10.1007/BF02480272}
#'
#'   Mohler, G. O., Short, M. B., Brantingham, P. J., Schoenberg, F. P. and
#'   Tita, G. E. (2011). Self-exciting point process modeling of crime.
#'   *Journal of the American Statistical Association*, 106(493), 100-108.
#'   \doi{10.1198/jasa.2011.ap09546}
#' @seealso [hwkres()] for the residuals, [hwksim()] to simulate.
#' @examples
#' x <- hwksim(c(a0 = log(0.5), eta = 0.5, beta = 2), 200, seed = 1)
#' fit <- hwkfit(x$t, 200)
#' fit
#' hwkfit(x$t, 200, kernel = c("exp", "lomax"))$cmp
#' @export
hwkfit <- function(
  times,
  hor,
  kernel = "exp",
  base = "const",
  per = 365.25,
  start = NULL,
  res = NULL,
  seed = NULL,
  maxit = 500
) {
  .Call(C_clab_r_hwkfit, times, hor, kernel, base, per, start, res, seed, maxit)
}

#' @export
print.hwkfit <- function(x, ...) {
  cat("hwkfit: ", x[["n"]], " events in [0, ", format(x[["hor"]]), "]; ",
      x[["kernel"]], " kernel, ", x[["base"]], " baseline\n", sep = "")
  print(data.frame(est = x[["par"]], se = x[["se"]]))
  b <- x[["branch"]]
  cat("branching ratio ", format(b[["eta"]]), ", mean cluster size ",
      format(b[["size"]]), "; log-likelihood ", format(x[["loglik"]]),
      ", AIC ", format(x[["aic"]]), "\n", sep = "")
  if (!x[["conv"]]) cat("the optimiser stopped before the gradient vanished\n")
  if (nrow(x[["cmp"]]) > 1L) print(x[["cmp"]], row.names = FALSE)
  invisible(x)
}

#' Hawkes process: time-rescaling residuals
#'
#' Transforms event times by the compensator of a Hawkes process,
#' \eqn{\Lambda(t) = \int_0^t \lambda}: under the model the gaps
#' \eqn{\tau_i = \Lambda(t_i) - \Lambda(t_{i-1})} (with \eqn{t_0 = 0}) are
#' independent standard exponentials, so `u = 1 - exp(-tau)` are uniform
#' on `[0, 1]`. The Kolmogorov-Smirnov test compares `u` with the uniform,
#' with R's exact p-value below 100 events without ties and the limiting
#' distribution otherwise, as [stats::ks.test()] does. `plot()` draws the
#' uniform quantile plot with the test's 95% band.
#'
#' @param times A fit from [hwkfit()], or event times as there.
#' @param hor,kernel,base,per As for [hwkfit()] (not read with a fit); one
#'   kernel and one baseline.
#' @param par The parameters, as `par` of [hwkfit()] (not read with a fit).
#' @return An object of class `hwkres`, a list: `tau`, the rescaled gaps;
#'   `u`, their uniform transforms; `d` and `p`, the Kolmogorov-Smirnov
#'   statistic and p-value; and `n`.
#' @references
#'   Ogata, Y. (1988). Statistical models for earthquake occurrences and
#'   residual analysis for point processes. *Journal of the American
#'   Statistical Association*, 83(401), 9-27.
#'   \doi{10.1080/01621459.1988.10478560}
#' @examples
#' x <- hwksim(c(a0 = log(0.5), eta = 0.5, beta = 2), 200, seed = 1)
#' r <- hwkres(hwkfit(x$t, 200))
#' r
#' plot(r)
#' @export
hwkres <- function(
  times,
  hor = NULL,
  par = NULL,
  kernel = "exp",
  base = "const",
  per = 365.25
) {
  if (inherits(times, "hwkfit")) {
    f <- times
    return(.Call(C_clab_r_hwkres, f[["times"]], f[["hor"]], f[["par"]],
                 f[["kernel"]], f[["base"]], f[["per"]]))
  }
  .Call(C_clab_r_hwkres, times, hor, par, kernel, base, per)
}

#' @export
print.hwkres <- function(x, ...) {
  cat("hwkres: ", x[["n"]], " rescaled gaps, mean ", format(mean(x[["tau"]])),
      " (1 under the model)\nKolmogorov-Smirnov D ", format(x[["d"]]),
      ", p ", format(x[["p"]]), "\n", sep = "")
  invisible(x)
}

#' @rdname hwkres
#' @param x An `hwkres` object.
#' @param ... Passed to [graphics::plot()].
#' @export
plot.hwkres <- function(x, ...) {
  n <- x[["n"]]
  q <- (seq_len(n) - 0.5) / n
  b <- 1.36 / sqrt(n)
  graphics::plot(q, sort(x[["u"]]), xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "uniform quantile", ylab = "rescaled time quantile",
                 ...)
  graphics::abline(0, 1)
  graphics::abline(b, 1, lty = 2)
  graphics::abline(-b, 1, lty = 2)
  invisible(x)
}

#' Hawkes process: simulation
#'
#' Simulates a Hawkes process on `[0, hor]` by its branching construction
#' (Hawkes and Oakes 1974): background events from the baseline (by
#' thinning for the `"cycle"` baseline), and from each event a Poisson
#' number, mean `eta`, of offspring at lags drawn from the kernel, and so on
#' over the generations. With `win` the process is in space and time:
#' background events fall uniformly in the rectangle, the baseline is then
#' per unit area, and each offspring lands at a normal offset (standard
#' deviation `sigma`) from its parent; offspring outside the window or past
#' `hor` are not recorded and have no offspring of their own.
#'
#' @param par A fit from [hwkfit()], or parameters as `par` of [hwkfit()].
#' @param hor The end of the window: one positive number.
#' @param kernel,base,per As for [hwkfit()]; one kernel and one baseline
#'   (not read with a fit).
#' @param win `NULL`, or the rectangle `c(xmin, xmax, ymin, ymax)` of a
#'   space-time process.
#' @param sigma The spread of offspring in space with `win`: one positive
#'   number.
#' @param seed `NULL` to draw the seed from R's generator, or one whole
#'   number from 0.
#' @param nmax The most events simulated: one whole number from 1 (more is
#'   an error).
#' @return A data frame in time order: `t`, with `x` and `y` in space, and
#'   `gen`, the generation (0 for a background event).
#' @references
#'   Hawkes, A. G. and Oakes, D. (1974). A cluster process representation of
#'   a self-exciting process. *Journal of Applied Probability*, 11(3),
#'   493-503. \doi{10.2307/3212693}
#' @examples
#' hwksim(c(a0 = log(0.2), eta = 0.6, beta = 1), 50, seed = 2)
#' hwksim(c(a0 = log(0.01), eta = 0.5, beta = 1), 50,
#'        win = c(0, 10, 0, 10), sigma = 0.5, seed = 2)
#' @export
hwksim <- function(
  par,
  hor,
  kernel = "exp",
  base = "const",
  per = 365.25,
  win = NULL,
  sigma = NULL,
  seed = NULL,
  nmax = 1e5
) {
  if (inherits(par, "hwkfit")) {
    f <- par
    return(.Call(C_clab_r_hwksim, f[["par"]], hor, f[["kernel"]], f[["base"]],
                 f[["per"]], win, sigma, seed, nmax))
  }
  .Call(C_clab_r_hwksim, par, hor, kernel, base, per, win, sigma, seed, nmax)
}
