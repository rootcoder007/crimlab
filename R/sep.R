#' Bayes outbreak detection
#'
#' The Bayes detector of `surveillance::algo.bayes()` (Riebler 2004;
#' Hoehle 2007) for a series of counts: at each time `t` the reference counts
#' are the `w` before it (`acty = TRUE`) and, in each of the `b` years before,
#' the `2 w + 1` counts around the same time of year (`freq` counts per
#' year); missing counts are left out. With `s` their sum and `n` their
#' number, the posterior predictive distribution of the count at `t` under
#' a Poisson model with a flat prior is negative binomial of size `s + 0.5`
#' and probability `n / (n + 1)`; its `1 - alpha` quantile is the upper
#' bound, and a count above it raises an alarm. The same series may be
#' calls for service, emergency visits or reported crimes of one kind.
#'
#' @param x The counts: numbers from 0, `NA` for a missing count.
#' @param freq The counts per year: one whole number from 1 (52 for weeks).
#' @param b The number of past years: one whole number from 0.
#' @param w The half-width of the window: one whole number from 0.
#' @param acty `TRUE` to use the `w` counts just before `t`.
#' @param alpha The level: one number strictly between 0 and 1.
#' @param at `NULL` for every time from `b freq + w + 1` on, or the times
#'   to evaluate (from 1).
#' @return A data frame with one row per time: `t`, `ub` (the upper bound)
#'   and `alarm` (`NA` where the count is missing).
#' @references Hoehle, M. (2007). surveillance: an R package for the
#'   monitoring of infectious diseases. *Computational Statistics*, 22(4),
#'   571-582. \doi{10.1007/s00180-007-0074-8}
#'
#'   Riebler, A. (2004). *Empirischer Vergleich von statistischen Methoden
#'   zur Ausbruchserkennung bei Surveillance Daten*. Bachelor's thesis,
#'   University of Munich.
#' @examples
#' bayout(c(3, 5, 2, 4, 6, 3, 4, 12), w = 6, at = 8)
#' @family spatial epidemiology
#' @export
bayout <- function(
  x,
  freq = 52,
  b = 0,
  w = 6,
  acty = TRUE,
  alpha = 0.05,
  at = NULL
) {
  .Call(C_clab_r_bayout, x, freq, b, w, acty, alpha, at)
}

#' Exposure within a buffer
#'
#' For each receptor (a home, a school), the number of sources (bars,
#' polluting sites, hot spots) within distance `r`, their summed weights
#' and the distance to the nearest: the buffer measure of exposure
#' assessment. Distances are planar, or great-circle metres when `lonlat =
#' TRUE` (the coordinates then longitude and latitude in degrees).
#'
#' @inheritParams cgstat
#' @param src The sources: a matrix or data frame of two coordinate
#'   columns, as `x`.
#' @param r The radius: one number from 0, in the unit of the coordinates
#'   (metres when `lonlat = TRUE`).
#' @param w `NULL` for a weight of 1 per source, or one weight from 0 per
#'   source.
#' @param lonlat `TRUE` when the receptors and sources are longitude and
#'   latitude in degrees.
#' @return A data frame with one row per receptor: `count`, `wsum` and
#'   `near`.
#' @references Nuckols, J. R., Ward, M. H. and Jarup, L. (2004). Using
#'   geographic information systems for exposure assessment in
#'   environmental epidemiology studies. *Environmental Health
#'   Perspectives*, 112(9), 1007-1015. \doi{10.1289/ehp.6738}
#' @examples
#' rec <- cbind(x = c(0, 1, 5), y = c(0, 1, 5))
#' src <- cbind(x = c(0.5, 4), y = c(0.5, 4))
#' bufexp(rec, src = src, r = 1.5)
#' @family spatial epidemiology
#' @export
bufexp <- function(x, y = NULL, src, r, w = NULL, lonlat = FALSE) {
  .Call(C_clab_r_bufexp, x, y, src, r, w, lonlat)
}

#' Exposure by kernel
#'
#' The kernel-weighted density of sources at each receptor, `sum_j w_j
#' K(d_ij / h) / h^2`, with the bivariate normal kernel `exp(-u^2 / 2) / (2
#' pi)` or the quartic kernel `3 / pi (1 - u^2)^2` (for `u < 1`, 0 beyond)
#' of Bailey and Gatrell (1995). Planar coordinates.
#'
#' @inheritParams bufexp
#' @param h The bandwidth: one positive number.
#' @param kern `"gaussian"` or `"quartic"`.
#' @return A numeric vector, one value per receptor.
#' @references Bailey, T. C. and Gatrell, A. C. (1995). *Interactive Spatial
#'   Data Analysis*. Longman, Harlow.
#' @examples
#' kerexp(rbind(c(0.2, 0.3)), src = rbind(c(0.25, 0.35), c(0.6, 0.8)),
#'        h = 0.5, kern = "quartic")
#' @family spatial epidemiology
#' @export
kerexp <- function(x, y = NULL, src, h, w = NULL, kern = "gaussian") {
  .Call(C_clab_r_kerexp, x, y, src, h, w, kern)
}

#' The Leroux precision matrix
#'
#' The precision of the area effects in the model of Leroux, Lei and
#' Breslow (2000), `tau (rho Q + (1 - rho) I)`, with `Q = D - W` the
#' Laplacian of the neighbour weights: independent effects at `rho = 0`,
#' the intrinsic conditional autoregression at `rho = 1`.
#'
#' @param w The neighbours: [spwt()]'s list, an spdep `listw` or an
#'   `n x n` matrix of weights from 0 (its non-zero entries off the
#'   diagonal).
#' @param rho The spatial parameter: one number from 0 to 1.
#' @param tau The precision: one positive number.
#' @return The `n x n` precision matrix.
#' @references Leroux, B. G., Lei, X. and Breslow, N. (2000). Estimation of
#'   disease rates in small areas: a new mixed model for spatial
#'   dependence. In M. E. Halloran and D. Berry (eds), *Statistical Models
#'   in Epidemiology, the Environment, and Clinical Trials*, 179-191.
#'   Springer, New York. \doi{10.1007/978-1-4612-1284-3_4}
#' @examples
#' lerprc(rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0)), 0.5, 2)
#' @family spatial epidemiology
#' @export
lerprc <- function(w, rho, tau = 1) {
  .Call(C_clab_r_lerprc, w, rho, tau)
}

#' The BYM2 area structure
#'
#' The scaled structure of the BYM2 model of Riebler et al. (2016): the
#' generalised inverse `G` of the Laplacian `Q` of the neighbour weights
#' (from its eigenvalues above `tol` times the largest), the scaling
#' factor `s`, the geometric mean of the diagonal of `G`, so that the
#' structured effect has a generalised variance of 1, and the covariance
#' of the area effects `((1 - phi) I + phi G / s) / tau`. The graph must be
#' connected and every area must have a neighbour.
#'
#' @inheritParams lerprc
#' @param phi The share of the variance that is structured: one number from
#'   0 to 1.
#' @param tol The relative eigenvalue tolerance: one number strictly
#'   between 0 and 1.
#' @return A list: `scale`, `ginv` and `cov`.
#' @references Riebler, A., Sorbye, S. H., Simpson, D. and Rue, H. (2016).
#'   An intuitive Bayesian spatial model for disease mapping that accounts
#'   for scaling. *Statistical Methods in Medical Research*, 25(4),
#'   1145-1165. \doi{10.1177/0962280216660421}
#' @examples
#' w <- rbind(c(0, 1, 1, 0), c(1, 0, 1, 0), c(1, 1, 0, 1), c(0, 0, 1, 0))
#' bym2sc(w, 2, 0.4)$scale
#' @family spatial epidemiology
#' @export
bym2sc <- function(w, tau = 1, phi = 0.5, tol = 1e-10) {
  .Call(C_clab_r_bym2sc, w, tau, phi, tol)
}

#' Ecological regression with expected counts
#'
#' Regression of area counts `y` on covariates with the offset `log(e)`, `e`
#' the expected counts (from indirect standardisation, say), so that the
#' coefficients act on the standardised ratio: Poisson by iteratively
#' reweighted least squares; negative binomial, alternating Newton steps on
#' the score of `theta` with reweighted fits of the coefficients; or
#' zero-inflated Poisson (Lambert 1992) by EM, the zero part logistic on an
#' intercept and the columns of `z`.
#'
#' @param y The counts: numbers from 0.
#' @param e The expected counts: positive numbers, one per count.
#' @param x `NULL`, or the covariates: a numeric matrix (or vector) with one
#'   row per count. An intercept is added, so a column that does not vary
#'   is refused.
#' @param family `"poisson"`, `"negbin"` or `"zip"`.
#' @param z For `"zip"`: `NULL`, or the covariates of the zero part, as `x`.
#' @param tol The convergence tolerance of `theta` and of the zero-inflated
#'   log-likelihood: one number strictly between 0 and 1.
#' @param maxit The rounds: one whole number from 1.
#' @return A list of class `ecolrg`: `coef` (the intercept first), `theta`
#'   (negative binomial), `zcoef` (zero-inflated), `fitted` (the means of
#'   the count part) and `loglik`.
#' @references Lambert, D. (1992). Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
#'   \doi{10.2307/1269547}
#'
#'   Lawson, A. B. (2018). *Bayesian Disease Mapping*, third edition.
#'   Chapman and Hall/CRC, Boca Raton. \doi{10.1201/9781351271769}
#' @examples
#' # simulated counts and expected counts
#' set.seed(1)
#' x <- rnorm(30)
#' e <- runif(30, 3, 8)
#' y <- rpois(30, e * exp(0.2 + 0.4 * x))
#' ecolrg(y, e, x)
#' @family spatial epidemiology
#' @export
ecolrg <- function(
  y,
  e,
  x = NULL,
  family = "poisson",
  z = NULL,
  tol = 1e-12,
  maxit = 500
) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.data.frame(z)) {
    z <- as.matrix(z)
  }
  .Call(C_clab_r_ecolrg, y, e, x, family, z, tol, maxit)
}

#' @export
print.ecolrg <- function(x, ...) {
  cat(
    "ecolrg: ",
    length(x[["fitted"]]),
    " counts, log-likelihood ",
    format(signif(x[["loglik"]], 6)),
    "\n",
    sep = ""
  )
  print(x[["coef"]])
  invisible(x)
}

#' The spatial encoding of graph attention
#'
#' Turns a matrix of hop distances between the nodes of a network (areas,
#' places, people) into the bias that Graphormer (Ying et al. 2021) adds to
#' the attention logits: a learnable value per distance, `b[d + 1]` (the
#' last for every larger distance), and `ub` for pairs that cannot reach
#' each other (distance `-1`). A bias inside the softmax keeps distant
#' nodes reachable but discouraged.
#'
#' @param d The hop distances: a matrix of whole numbers from 0, `-1` for
#'   an unreachable pair.
#' @param b The bias table, indexed by distance from 0.
#' @param ub The bias of an unreachable pair: one number.
#' @return A list: `bias` (a matrix the shape of `d`) and `ub`.
#' @references Ying, C., Cai, T., Luo, S., Zheng, S., Ke, G., He, D.,
#'   Shen, Y. and Liu, T.-Y. (2021). Do transformers really perform badly
#'   for graph representation? *Advances in Neural Information Processing
#'   Systems*, 34, 28877-28888. \doi{10.48550/arXiv.2106.05234}
#' @examples
#' d <- rbind(c(0, 1, 2), c(1, 0, -1), c(2, -1, 0))
#' spbias(d, b = c(1, 0.5, 0.25))
#' @family spatial epidemiology
#' @export
spbias <- function(d, b, ub = -10) {
  .Call(C_clab_r_spbias, d, b, ub)
}

#' The metapopulation SIR model
#'
#' The deterministic SIR model on patches (Keeling and Rohani 2008, sec.
#' 7.3): the force of infection on patch `i` is `beta sum_j C_ij I_j /
#' N_j` with the coupling `C = (1 - kappa) I + kappa W~`, `W~` the
#' row-standardised neighbour weights, so that a share `kappa` of the
#' pressure comes from neighbours; recovery is at rate `gamma`. Integrated
#' by the classical fourth-order Runge-Kutta method in steps of `dt`. The
#' same equations describe the spread of an outbreak of violence between
#' neighbourhoods in contagion models.
#'
#' @param i0 The infected at time 0, one number from 0 per patch.
#' @param w `NULL` for isolated patches, or the neighbours, as in
#'   [lerprc()].
#' @param pop `NULL` for 1000 per patch, or the populations (positive, at
#'   least `i0`).
#' @param beta,gamma The transmission and recovery rates: numbers from 0.
#' @param kappa The share of the pressure from neighbours: from 0 to 1.
#' @param tmax,dt The horizon and the step: positive; at most 1e6 steps.
#' @return A list of class `spsir`: `attack` (the final attack rate, `sum
#'   R / sum N`), `r0` (`beta / gamma`), `t` (the times) and `S`, `I`, `R`
#'   (matrices: a row per time, a column per patch).
#' @references Keeling, M. J. and Rohani, P. (2008). *Modeling Infectious
#'   Diseases in Humans and Animals*. Princeton University Press,
#'   Princeton. \doi{10.1515/9781400841035}
#' @examples
#' w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
#' spsir(c(10, 0, 0), w, rep(1000, 3), beta = 0.6, gamma = 0.2,
#'       tmax = 60, dt = 0.5)
#' @family spatial epidemiology
#' @export
spsir <- function(
  i0,
  w = NULL,
  pop = NULL,
  beta = 0.5,
  gamma = 0.2,
  kappa = 0.1,
  tmax = 100,
  dt = 0.1
) {
  .Call(C_clab_r_spsir, i0, w, pop, beta, gamma, kappa, tmax, dt)
}

#' @export
print.spsir <- function(x, ...) {
  cat(
    "spsir: ",
    ncol(x[["S"]]),
    " patches to time ",
    format(x[["t"]][length(x[["t"]])]),
    ", attack rate ",
    format(signif(x[["attack"]], 4)),
    ", R0 ",
    format(signif(x[["r0"]], 4)),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' Spatial frailty and cure models
#'
#' `sfrail()` fits the Weibull proportional hazards model with log-frailties
#' of the areas under an intrinsic conditional autoregression,
#' `h(t) = rho t^(rho - 1) exp(x` `beta + u_area)` (Banerjee, Wall and Carlin
#' 2003);
#' `sfcure()` the promotion-time cure model (Chen, Ibrahim and Sinha 1999),
#' `S(t) = exp(-theta (1 -` `exp(-G(t))))`, `theta = exp(x beta + u_area)` and
#' `G(t) = exp(c + rho log t)`, whose cure fraction is `exp(-theta)`. Both
#' maximise the penalised likelihood, the penalty `tau u'Qu / 2` with `Q`
#' the Laplacian of the neighbour weights, by Newton-Raphson from 0 on the
#' analytic gradient with a Hessian by central differences, the step halved
#' while the objective falls, until a step is below `tol` (200 iterations
#' at most). The standard errors come from the last Hessian; the first is
#' that of the log shape.
#'
#' @param time The times: positive numbers.
#' @param event The events: 1 (or `TRUE`) for a failure, 0 for a censored
#'   time.
#' @param x `NULL`, or the covariates (no intercept: the area effects hold
#'   the level), a numeric matrix (or vector) with one row per subject.
#' @param reg The area of each subject: whole numbers from 1 to the number
#'   of areas of `w`.
#' @inheritParams lerprc
#' @param tau The precision of the area effects: one number from 0.
#' @param tol The convergence tolerance: one number strictly between 0 and
#'   1.
#' @return For `sfrail()`, a list of class `sfrail`: `pll` (the penalised
#'   log-likelihood), `shape`, `coef`, `frail`, `se`, `iter`, `tau`. For
#'   `sfcure()`, of class `sfcure`: `pll`, `shape`, `lscale`, `coef`,
#'   `effect`, `se`, `cure` (each subject's cure fraction), `iter`, `tau`.
#'   [summary()] gives the table of estimates.
#' @references Banerjee, S., Wall, M. M. and Carlin, B. P. (2003). Frailty
#'   modeling for spatially correlated survival data, with application to
#'   infant mortality in Minnesota. *Biostatistics*, 4(1), 123-142.
#'   \doi{10.1093/biostatistics/4.1.123}
#'
#'   Chen, M.-H., Ibrahim, J. G. and Sinha, D. (1999). A new Bayesian model
#'   for survival data with a surviving fraction. *Journal of the American
#'   Statistical Association*, 94(447), 909-919.
#'   \doi{10.1080/01621459.1999.10474196}
#' @examples
#' # simulated survival times in three areas
#' w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
#' i <- 0:29
#' x <- sin(i)
#' t <- 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * x - 0.3 * (i %% 3))
#' f <- sfrail(t, as.numeric(i %% 5 != 0), x, i %% 3 + 1, w)
#' f
#' summary(f)
#' @family spatial epidemiology
#' @export
sfrail <- function(time, event, x = NULL, reg, w, tau = 1, tol = 1e-10) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  .Call(C_clab_r_sfrail, time, event, x, reg, w, tau, tol)
}

#' @rdname sfrail
#' @export
sfcure <- function(time, event, x = NULL, reg, w, tau = 1, tol = 1e-10) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  .Call(C_clab_r_sfcure, time, event, x, reg, w, tau, tol)
}

#' @export
print.sfrail <- function(x, ...) {
  cat(
    "sfrail: Weibull frailty model, shape ",
    format(signif(x[["shape"]], 4)),
    ", penalised log-likelihood ",
    format(signif(x[["pll"]], 6)),
    "\n",
    sep = ""
  )
  print(c(x[["coef"]], x[["frail"]]))
  invisible(x)
}

#' @export
summary.sfrail <- function(object, ...) {
  b <- c(log(object[["shape"]]), object[["coef"]], object[["frail"]])
  data.frame(
    term = .sfnames(object, "logshape", "frail"),
    est = b,
    se = object[["se"]],
    z = b / object[["se"]]
  )
}

#' @export
print.sfcure <- function(x, ...) {
  cat(
    "sfcure: promotion-time cure model, shape ",
    format(signif(x[["shape"]], 4)),
    ", penalised log-likelihood ",
    format(signif(x[["pll"]], 6)),
    "\n",
    sep = ""
  )
  print(summary(x[["cure"]]))
  invisible(x)
}

#' @export
summary.sfcure <- function(object, ...) {
  b <- c(
    log(object[["shape"]]),
    object[["lscale"]],
    object[["coef"]],
    object[["effect"]]
  )
  data.frame(
    term = c("logshape", .sfnames(object, "lscale", "effect")),
    est = b,
    se = object[["se"]],
    z = b / object[["se"]]
  )
}

#' Eigenvector spatial filtering
#'
#' Adds Moran eigenvectors of the neighbour weights to the regression of
#' `v` on an intercept and `x` by forward selection, so that the residuals
#' lose their spatial autocorrelation (Griffith 2003; Tiefelsdorf and
#' Griffith 2007). The Moran eigenvectors are the eigenvectors of the
#' symmetrised weights `(W + W') / 2` orthogonal to the constant (computed
#' in the Helmert basis), each with Moran's I `n / S0 lambda`; the
#' candidates are those with Moran's I above `-1 / (n - 1)` (`"positive"`,
#' the usual choice), below it, or all. Selection by `method`: `"moran"`
#' takes the vector that most reduces the absolute residual Moran's I until
#' it is below `tol`; `"aic"` and `"press"` (the leave-one-out prediction
#' error) while they fall; `"r2"` while R-squared rises by `tol` or more.
#' At least two residual degrees of freedom are kept.
#'
#' @param v The outcome: finite numbers, one per area, not all the same.
#' @param w The neighbours, as in [lerprc()], for as many areas as `v`.
#' @param x `NULL`, or covariates: a numeric matrix (or vector) with one
#'   row per area. An intercept is added, so a column that does not vary
#'   is refused.
#' @param method `"moran"`, `"aic"`, `"press"` or `"r2"`.
#' @param cand `"positive"`, `"negative"` or `"all"`.
#' @param tol The stopping tolerance: one positive number.
#' @param maxv `NULL`, or the largest number of vectors.
#' @return A list of class `evfilt`: `sel` (the selected eigenvectors in
#'   the order they were added, numbered by decreasing eigenvalue), `coef`
#'   (the intercept, `x`, then the vectors), `imoran`
#'   (the residual Moran's I), `aic`, `r2`, `fr2` (the gain in R-squared
#'   from the filter), `resid` and `mi` (the Moran's I of every
#'   eigenvector).
#' @references Griffith, D. A. (2003). *Spatial Autocorrelation and Spatial
#'   Filtering*. Springer, Berlin. \doi{10.1007/978-3-540-24806-4}
#'
#'   Tiefelsdorf, M. and Griffith, D. A. (2007). Semiparametric filtering of
#'   spatial autocorrelation: the eigenvector approach. *Environment and
#'   Planning A*, 39(5), 1193-1221. \doi{10.1068/a37378}
#' @examples
#' w <- 1 * (abs(outer(1:8, 1:8, "-")) == 1)
#' v <- c(1, 2, 2.5, 4, 3.5, 3, 1.5, 1)
#' evfilt(v, w, method = "aic")
#' @family spatial epidemiology
#' @export
evfilt <- function(
  v,
  w,
  x = NULL,
  method = "moran",
  cand = "positive",
  tol = 0.1,
  maxv = NULL
) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  .Call(C_clab_r_evfilt, v, w, x, method, cand, tol, maxv)
}

#' @export
print.evfilt <- function(x, ...) {
  cat(
    "evfilt: ",
    length(x[["sel"]]),
    " eigenvectors selected, residual ",
    "Moran's I ",
    format(signif(x[["imoran"]], 4)),
    ", R-squared ",
    format(signif(x[["r2"]], 4)),
    "\n",
    sep = ""
  )
  invisible(x)
}

.sfnames <- function(x, a, b) {
  c(
    a,
    if (length(x[["coef"]])) paste0("x", seq_along(x[["coef"]])),
    paste0(b, seq_along(x[[b]]))
  )
}
