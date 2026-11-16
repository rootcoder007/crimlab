# Counts over time and across categories: thin wrappers over the kernels
# of src/tsd.cpp (src/svk.h).

#' ARIMA(1,1,1) forecast of a count series
#'
#' Fits ARIMA(1,1,1) to a series of counts (monthly incidents, calls for
#' service, admissions) by exact maximum likelihood, the Kalman filter on
#' the state-space form with the conditional-sum-of-squares fit as its
#' start, as [stats::arima()] does by default, and forecasts `h` periods
#' ahead with their standard errors. AIC is \eqn{-2 \log L + 2(k + 1)};
#' BIC is \eqn{\mathrm{AIC} + (k + 1)(\log n - 2)}, \eqn{k} the ARMA
#' coefficients.
#'
#' @param y The counts, in time order, 24 at least.
#' @param h The periods to forecast, at most 67108864 (2^26).
#' @return A list of class `arimfc`: `fc` and `se` (forecasts and their
#'   standard errors), `coef` (ar1, ma1), `sigma2`, `aic`, `bic` and `n`.
#' @references Box, G. E. P. and Jenkins, G. M. (1970). *Time Series
#'   Analysis: Forecasting and Control*. Holden-Day.
#'
#'   Gardner, G., Harvey, A. C. and Phillips, G. D. A. (1980). Algorithm
#'   AS 154: an algorithm for exact maximum likelihood estimation of
#'   autoregressive-moving average models by means of Kalman filtering.
#'   *Journal of the Royal Statistical Society C*, 29(3), 311-322.
#' @examples
#' # 48 simulated months of counts
#' set.seed(1)
#' y <- rpois(48, 50 + 10 * sin(1:48 / 6))
#' arimfc(y, h = 6)$fc
#' @family counts over time
#' @export
arimfc <- function(y, h = 12) .sv("arimfc", list(y = y, h = h), "tsfit")

#' Seasonal ARIMA forecast with a hold-out check
#'
#' Fits SARIMA\eqn{(p,d,q)(P,D,Q)_s} (by the method of [arimfc()]) to
#' all but the last `h` values of a count series, forecasts those `h`,
#' and scores the forecasts against them: the mean absolute percentage
#' error (each error over the larger of the count and 1) and the root
#' mean square error. The model must difference the series once at
#' least (`order[2]` or `seasonal[2]` above 0), so it has no mean.
#'
#' @param y The counts, in time order, 36 at least.
#' @param h The values held out.
#' @param order \eqn{(p, d, q)}.
#' @param seasonal \eqn{(P, D, Q, s)}, \eqn{s} the period (12 for
#'   months).
#' @return A list of class `sarfc`: `fc`, `actual`, `coef`, `aic`, `bic`
#'   (with the training length), `mape` (percent), `rmse` and `ntrain`.
#' @references Box, G. E. P. and Jenkins, G. M. (1970). *Time Series
#'   Analysis: Forecasting and Control*. Holden-Day.
#'
#'   Gardner, G., Harvey, A. C. and Phillips, G. D. A. (1980). Algorithm
#'   AS 154: an algorithm for exact maximum likelihood estimation of
#'   autoregressive-moving average models by means of Kalman filtering.
#'   *Journal of the Royal Statistical Society C*, 29(3), 311-322.
#' @examples
#' set.seed(1)
#' y <- rpois(60, 50 + 10 * sin(2 * pi * (1:60) / 12))
#' sarfc(y, h = 6)[c("mape", "rmse")]
#' @family counts over time
#' @export
sarfc <- function(y, h = 12, order = c(1, 1, 1), seasonal = c(0, 1, 1, 12)) {
  .sv("sarfc", list(y = y, h = h, order = order, seasonal = seasonal),
      "tsfit")
}

#' Correlation of counts per unit between categories
#'
#' Counts the events of each category (offence type, call type) by unit
#' (neighbourhood, beat, district), over the union of the units with a 0
#' where a category has none, and gives the Pearson correlations of those
#' counts between every pair of categories: whether the same places score
#' high on both.
#'
#' @param units A named list, one element per category: the unit of each
#'   event (`NA` dropped).
#' @return A list of class `catcor`: `r` (categories by categories, `NA`
#'   for a category of one count everywhere), `counts` (units by
#'   categories), `units` and `cats` (sorted).
#' @references Pearson, K. (1896). Mathematical contributions to the
#'   theory of evolution. III. Regression, heredity, and panmixia.
#'   *Philosophical Transactions of the Royal Society A*, 187, 253-318.
#'
#' @examples
#' # simulated events of three categories in 20 areas
#' set.seed(1)
#' a <- sprintf("%02d", 1:20)
#' u <- list(theft = sample(a, 200, TRUE), assault = sample(a, 80, TRUE),
#'           fraud = sample(a, 50, TRUE))
#' catcor(u)$r
#' @family counts over time
#' @export
catcor <- function(units) {
  if (!is.list(units)) {
    stop("catcor: `units` must be a named list of categories.", call. = FALSE)
  }
  u <- lapply(units, function(e) as.character(e[!is.na(e)]))
  r <- .sv("catcor", list(units = u), "tsfit")
  dimnames(r$r) <- list(r$cats, r$cats)
  dimnames(r$counts) <- list(r$units, r$cats)
  r
}

#' Ornstein-Uhlenbeck (Langevin) fit and simulation of a count series
#'
#' Treats a daily count series as a noisy mean-reverting intensity,
#' \eqn{dX_t = \theta (\mu - X_t) dt + \sigma dW_t}: \eqn{\theta},
#' \eqn{\mu} and \eqn{\sigma} from the least-squares regression of the
#' first differences on the lagged values (\eqn{\theta = -b}, at least
#' \eqn{10^{-6}}, \eqn{\mu = -a/b}, \eqn{\sigma} the residual standard
#' deviation), then `npath` Euler-Maruyama paths of `days / dt` steps from
#' the last value. Draws come from `seed`, as [set.seed()] takes it.
#'
#' @param y The counts, in time order, three at least.
#' @param npath The paths.
#' @param days The horizon, in the time unit of `y`.
#' @param dt The step.
#' @inheritParams bamscl
#' @return A list of class `ousim`: `theta`, `mu`, `sigma`, `halflife`
#'   (\eqn{\log 2 / \theta}), `q` (the 5%, 50% and 95% points of the last
#'   step), `paths` (paths by steps) and `n`.
#' @references Uhlenbeck, G. E. and Ornstein, L. S. (1930). On the theory
#'   of the Brownian motion. *Physical Review*, 36(5), 823-841.
#'   \doi{10.1103/PhysRev.36.823}
#'
#'   Kloeden, P. E. and Platen, E. (1992). *Numerical Solution of
#'   Stochastic Differential Equations*. Springer.
#'   \doi{10.1007/978-3-662-12616-5}
#' @examples
#' set.seed(1)
#' y <- rpois(200, 20)
#' ousim(y, npath = 20, days = 30)$q
#' @family counts over time
#' @export
ousim <- function(y, npath = 100, days = 365, dt = 1, seed = 42) {
  .sv("ousim", list(y = y, npath = npath, days = days, dt = dt, seed = seed),
      "tsfit")
}

#' Fokker-Planck evolution of a count density
#'
#' The density of the Ornstein-Uhlenbeck process fitted as in [ousim()],
#' evolved by the Fokker-Planck equation \eqn{\partial_t p = -\partial_x
#' (\theta (\mu - x) p) + \frac{1}{2} \sigma^2 \partial_x^2 p}, here with
#' the drift taken outside the derivative as in the explicit scheme
#' used: central differences on `ngrid` points from 0 to 1.5 times the
#' largest count plus 1, `nstep` steps of 0.05, edges copied from their
#' neighbours, the density kept non-negative and of total 1. It starts as
#' a normal curve at the last count with standard deviation
#' \eqn{\sigma + 1}. The stationary variance is \eqn{\sigma^2 / (2
#' \theta)}.
#'
#' @inheritParams ousim
#' @param ngrid The grid points, at most 67108864 (2^26).
#' @param nstep The steps.
#' @return A list of class `fpdens`: `theta`, `mu`, `sigma`, `grid`,
#'   `density`, `peak` (the grid point of the highest density) and
#'   `statvar`.
#' @references Risken, H. (1989). *The Fokker-Planck Equation*, 2nd ed.
#'   Springer.
#' @examples
#' set.seed(1)
#' f <- fpdens(rpois(200, 20))
#' f$peak
#' @family counts over time
#' @export
fpdens <- function(y, ngrid = 64, nstep = 200) {
  .sv("fpdens", list(y = y, ngrid = ngrid, nstep = nstep), "tsfit")
}

#' Lotka-Volterra rates of a yearly count series
#'
#' Treats yearly counts as the prey \eqn{x(t)} and their three-year
#' trailing mean as a stand-in predator \eqn{y(t)} (an enforcement series
#' can take its place) in \eqn{\dot x = \alpha x - \beta x y}, \eqn{\dot
#' y = \delta x y - \gamma y}: \eqn{\alpha} is the median log growth of
#' \eqn{x} plus 0.05, \eqn{\gamma} minus the median log growth of
#' \eqn{y} plus 0.05, \eqn{\beta = \alpha / \mathrm{median}(y)},
#' \eqn{\delta = \gamma / \mathrm{median}(x)}, and the linearised cycle
#' has period \eqn{2\pi / \sqrt{\alpha\gamma}} (D'Orsogna and Perc 2015,
#' section 3.4). The rates describe the series; they are not a fitted
#' dynamical model.
#'
#' @param y The yearly counts (counts below 1 count as 1), five at least.
#' @param year `NULL` or the year of each count.
#' @return A list of class `lvfit`: `alpha`, `beta`, `gamma`, `delta`,
#'   `period`, `year`, `prey` and `pred`.
#' @references D'Orsogna, M. R. and Perc, M. (2015). Statistical physics
#'   of crime: a review. *Physics of Life Reviews*, 12, 1-21.
#'   \doi{10.1016/j.plrev.2014.11.001}
#' @examples
#' lvfit(c(410, 455, 430, 470, 520, 498, 515), year = 2017:2023)$period
#' @family counts over time
#' @export
lvfit <- function(y, year = NULL) {
  .sv("lvfit", list(y = y, year = year), "tsfit")
}
