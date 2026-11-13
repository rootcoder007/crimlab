#' ETAS: a self-exciting process with marks and the Omori law
#'
#' Fits the temporal epidemic-type aftershock sequence (ETAS) model of
#' Ogata (1988) to event times with marks (a magnitude, the number of
#' victims, the severity of an incident): \deqn{\lambda(t) = \mu + K
#' \sum_{t_j < t} e^{\alpha (m_j - m_0)} (t - t_j + c)^{-p},} so bigger
#' events trigger more, and the triggering decays as a power of the time
#' since (the modified Omori law). Without marks, or with constant ones,
#' \eqn{\alpha} is fixed at 0 and the model is a Hawkes process with a
#' Lomax kernel (that of [hwkfit()] with shape \eqn{p - 1} and scale
#' \eqn{c}). The log-likelihood over `[0, hor]` is maximised with its
#' analytic gradient by BFGS over the logarithms of \eqn{\mu}, \eqn{K},
#' \eqn{c} and \eqn{p - 1} (and \eqn{\alpha}), from two starts or `start`;
#' standard errors from the numerical Hessian of the gradient. `branch` is
#' the branching ratio \eqn{K \bar{e^{\alpha (m - m_0)}} c^{1 - p} / (p -
#' 1)}, the mean number of events each event triggers. On data with
#' little power-law tail the likelihood rises towards the exponential
#' limit of the kernel (`c` and `p` growing together): the fit then stops
#' with `conv` `FALSE`, and [hwkfit()]'s exponential kernel is the better
#' model.
#'
#' @param times Event times: numbers sorted in `[0, hor]`, at least 10.
#' @inheritParams hwksim
#' @param mag `NULL` or the mark of each event: finite numbers.
#' @param m0 `NULL` (the least mark) or the reference mark: one number.
#' @param start `NULL` or the starting `mu`, `k`, `alpha`, `c` and `p`.
#' @inheritParams hwkfit
#' @return An object of class `etasf`, a list: `par`, `se` and `cov` of
#'   (`mu`, `k`, `alpha`, `c`, `p`; `alpha`'s are `NA` when it is fixed);
#'   `loglik`, `aic`; `branch`; `conv`, `iter`; `afree`, whether `alpha`
#'   was estimated; `n` and `hor`.
#' @references
#'   Ogata, Y. (1988). Statistical models for earthquake occurrences and
#'   residual analysis for point processes. *Journal of the American
#'   Statistical Association*, 83(401), 9-27.
#'   \doi{10.1080/01621459.1988.10478560}
#' @seealso [hwkfit()], [mvhwk()], [sthwk()].
#' @examples
#' x <- hwksim(c(log(0.3), 0.5, 0.6, 1), 800, kernel = "lomax", seed = 5)
#' etasf(x$t, 800)
#' @export
etasf <- function(times, hor, mag = NULL, m0 = NULL, start = NULL,
                  maxit = 500) {
  .Call(C_clab_r_etasf, times, hor, mag, m0, start, maxit)
}

#' @export
print.etasf <- function(x, ...) {
  cat("etasf: ", x[["n"]], " events in [0, ", format(x[["hor"]]), "]",
      if (!x[["afree"]]) "; no marks, alpha fixed at 0", "\n", sep = "")
  print(data.frame(est = x[["par"]], se = x[["se"]]))
  cat("branching ratio ", format(x[["branch"]]), "; log-likelihood ",
      format(x[["loglik"]]), ", AIC ", format(x[["aic"]]), "\n", sep = "")
  invisible(x)
}

#' Multivariate Hawkes process: types that excite each other
#'
#' Fits a multivariate Hawkes process with exponential kernels to event
#' times of several types (offence types, gangs, areas): \deqn{\lambda_c(t)
#' = \mu_c + \beta \sum_l A_{cl} \sum_{t_j < t, \, type_j = l}
#' e^{-\beta (t - t_j)},} where \eqn{A_{cl}} is the expected number of
#' events of type `c` triggered by one event of type `l` (retaliation
#' between gangs, one crime type leading to another). The rates, the
#' matrix and the shared decay `beta` (or a given one) are fitted by
#' maximum likelihood with the analytic gradient (Ozaki's recursion, O(n
#' k)) by BFGS over their logarithms. The process is stationary when the
#' spectral radius `rad` of \eqn{A} is below 1; with one type it is the
#' exponential process of [hwkfit()] with `eta` \eqn{= A}.
#'
#' @param times Event times: numbers sorted in `[0, hor]`, at least 10.
#' @param type The type of each event: labels, at most 10 types.
#' @inheritParams hwksim
#' @param beta `NULL` to estimate the decay, or a given one: one positive
#'   number.
#' @param start `NULL` or the starting parameters in the order of `par`
#'   (positive).
#' @param maxit The most BFGS steps: one whole number from 1.
#' @return An object of class `mvhwk`, a list: `par`, `se` and `cov` (named
#'   `mu.<type>`, `a.<to>.<from>`, `beta`); `mu`, the rates; `a`, the
#'   matrix (rows the triggered type, columns the triggering type);
#'   `beta`; `rad`; `loglik`, `aic`; `conv`, `iter`; `n` and `hor`.
#' @references
#'   Hawkes, A. G. (1971). Spectra of some self-exciting and mutually
#'   exciting point processes. *Biometrika*, 58(1), 83-90.
#'   \doi{10.1093/biomet/58.1.83}
#'
#'   Linderman, S. W. and Adams, R. P. (2014). Discovering latent network
#'   structure in point process data. *Proceedings of the 31st
#'   International Conference on Machine Learning*, 32(2), 1413-1421.
#' @examples
#' a <- hwksim(c(log(0.3), 0.3, 1), 400, seed = 1)
#' b <- hwksim(c(log(0.2), 0.4, 1), 400, seed = 2)
#' d <- rbind(data.frame(t = a$t, type = "a"), data.frame(t = b$t, type = "b"))
#' d <- d[order(d$t), ]
#' mvhwk(d$t, d$type, 400)
#' @export
mvhwk <- function(times, type, hor, beta = NULL, start = NULL, maxit = 500) {
  .Call(C_clab_r_mvhwk, times, type, hor, beta, start, maxit)
}

#' @export
print.mvhwk <- function(x, ...) {
  cat("mvhwk: ", x[["n"]], " events of ", length(x[["mu"]]),
      " types; beta ", format(x[["beta"]]), ", spectral radius ",
      format(x[["rad"]]), "\nbackground rates:\n", sep = "")
  print(x[["mu"]])
  cat("excitation (row: triggered, column: triggering):\n")
  print(x[["a"]])
  invisible(x)
}

#' Space-time Hawkes process
#'
#' Fits the self-exciting point process of crime in space and time of
#' Mohler et al. (2011): \deqn{\lambda(t, x, y) = \mu + \eta \sum_{t_j <
#' t} \beta e^{-\beta (t - t_j)} \phi_\sigma(x - x_j, y - y_j),} a
#' constant background \eqn{\mu} per unit area and time over the window
#' and, after each event, `eta` expected offspring at exponential lags
#' (rate `beta`) and normal offsets (standard deviation `sigma`). The
#' log-likelihood takes the spatial kernel's mass as 1 inside the window
#' (events near its edge lose a little offspring mass outside); it is
#' maximised with its analytic gradient by BFGS over the logarithms of
#' `mu`, `beta`, `sigma` and the logit of `eta`, from two starts or
#' `start`. With `par` nothing is fitted: the log-likelihood is evaluated
#' there. `at` gives the conditional intensity at query points, from the
#' events strictly before each. [hwksim()] simulates from a fit.
#'
#' @param times Event times: numbers sorted in `[0, hor]`, at least 10.
#' @param x,y The places of the events: finite numbers in planar units.
#' @param hor The end of the window in time: one positive number.
#' @param win `NULL` (the events' bounding box) or the rectangle
#'   `c(xmin, xmax, ymin, ymax)` holding every event.
#' @param par `NULL` to fit, or the parameters `mu`, `eta`, `beta`,
#'   `sigma` at which to evaluate.
#' @param at `NULL` or query points: a data frame or list with `t`, `x`
#'   and `y`.
#' @param start `NULL` or the starting parameters, as `par`.
#' @inheritParams hwkfit
#' @return An object of class `sthwk`, a list: `par`, `se` and `cov`
#'   (`mu`, `eta`, `beta`, `sigma`); `loglik`, `aic`; `conv`, `iter`;
#'   `fit` (`FALSE` when `par` was given); `n`, `hor`, `win`, `area`; and
#'   `lam`, the intensities at `at` (`NULL` without it).
#' @references
#'   Mohler, G. O., Short, M. B., Brantingham, P. J., Schoenberg, F. P. and
#'   Tita, G. E. (2011). Self-exciting point process modeling of crime.
#'   *Journal of the American Statistical Association*, 106(493), 100-108.
#'   \doi{10.1198/jasa.2011.ap09546}
#'
#'   Reinhart, A. (2018). A review of self-exciting spatio-temporal point
#'   processes and their applications. *Statistical Science*, 33(3),
#'   299-318. \doi{10.1214/17-STS629}
#' @examples
#' s <- hwksim(c(log(0.02), 0.5, 1), 100, win = c(0, 10, 0, 10),
#'             sigma = 0.3, seed = 3)
#' f <- sthwk(s$t, s$x, s$y, 100, win = c(0, 10, 0, 10))
#' f
#' sthwk(s$t, s$x, s$y, 100, win = c(0, 10, 0, 10), par = f$par,
#'       at = data.frame(t = 50, x = 5, y = 5))$lam
#' @export
sthwk <- function(
  times,
  x,
  y,
  hor,
  win = NULL,
  par = NULL,
  at = NULL,
  start = NULL,
  maxit = 500
) {
  .Call(C_clab_r_sthwk, times, x, y, hor, win, par, at, start, maxit)
}

#' @export
print.sthwk <- function(x, ...) {
  cat("sthwk: ", x[["n"]], " events in [0, ", format(x[["hor"]]),
      "] over an area of ", format(x[["area"]]),
      if (!x[["fit"]]) "; evaluated at the given parameters", "\n",
      sep = "")
  print(data.frame(est = x[["par"]], se = x[["se"]]))
  cat("log-likelihood ", format(x[["loglik"]]), "\n", sep = "")
  invisible(x)
}
