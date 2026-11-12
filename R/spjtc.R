#' Mantel's space-time test
#'
#' The Pearson correlation, over all pairs of events, between their
#' spatial dissimilarity (the distance, or `1 / (distance + cs)` when
#' `cs > 0`, Mantel's reciprocal transform) and their temporal one
#' (`|t_i - t_j|`, or `1 / (|t_i - t_j| + ct)`), with a permutation test
#' that permutes the times over the fixed locations (Mantel 1967). Drawing
#' distances and intervals independently instead, as some descriptions do,
#' does not keep the dependence among the pairs of one distance matrix and
#' gives no test of stated size. Memory grows with the number of events,
#' not pairs; time with pairs times permutations.
#'
#' @inheritParams cgstat
#' @param t The time of each event: finite numbers (for dates, `as.numeric()`
#'   of them, in days or seconds).
#' @param cs,ct The constants of the reciprocal transforms: numbers from 0
#'   (0 uses the raw distance or interval).
#' @param nsim The number of permutations: one whole number from 0.
#' @inheritParams nnidx
#' @return A named numeric vector: `r` (the correlation), `pairs`, and the
#'   permutation p-values `pup` (`(1 + #{r_perm >= r}) / (nsim + 1)`, large
#'   r: events close in space are close in time when both are raw
#'   distances), `plo` (`r_perm <= r`) and `p2` (`|r_perm| >= |r|`); `NaN`
#'   without permutations. 3 to 46340 events.
#' @references Mantel, N. (1967). The detection of disease clustering and a
#'   generalized regression approach. *Cancer Research*, 27(2), 209-220.
#' @examples
#' set.seed(1)
#' x <- runif(60)
#' y <- runif(60)
#' mantst(x, y, t = x + rnorm(60, sd = 0.1), nsim = 199, seed = 1)
#' @export
mantst <- function(x, y = NULL, t, cs = 0, ct = 0, nsim = 999, seed = NULL,
                   lonlat = FALSE) {
  .Call(C_clab_r_mantst, x, y, t, cs, ct, nsim, seed, lonlat)
}

#' Spatial-temporal moving average
#'
#' For each event, the mean centre of the `span` events centred on it in
#' time order (ties in input order), the window shortened at the ends of
#' the series. It traces how the centre of a series of incidents moves
#' over time. It is a descriptive device of the crime-mapping literature
#' with no original paper to cite; the mean centre is that of [cgstat()].
#'
#' @inheritParams mantst
#' @param span The window: one odd whole number from 1.
#' @return A data frame with one row per event, in input order: `mx`, `my`
#'   (the window's mean centre, in degrees when `lonlat = TRUE`) and `n`
#'   (events in the window).
#' @examples
#' stmavg(x = 1:7, y = c(1, 1, 2, 3, 5, 8, 13), t = 7:1, span = 3)
#' @export
stmavg <- function(x, y = NULL, t, span = 5, lonlat = FALSE) {
  .Call(C_clab_r_stmavg, x, y, t, span, lonlat)
}

#' Distance-decay functions for journey to crime
#'
#' The decay of the chance of offending with distance from an offender's
#' base, in the forms of Brantingham and Brantingham (1981), Rossmo (2000)
#' and Levine (2009): `"linear"`, `f = max(0, A + B d)`, `par = c(A, B)`;
#' `"negexp"`, `f = A exp(-b d)`, `par = c(A, b)` with `b >= 0`;
#' `"normal"`, `f = A dnorm(d, mu, s)`, `par = c(A, mu, s)`; `"lognormal"`,
#' `f = A dlnorm(d, mu, s)`, `par = c(A, mu, s)` (the lognormal density of
#' `d`; some texts use `log d^2` and `1 / d^2`, which is not a density);
#' `"trunc"`, the truncated negative exponential, `f = P d / dp` up to the
#' peak `dp` and `P exp(-C (d - dp))` beyond, `par = c(P, dp, C)`; and
#' `"table"`, values `tab` at distances `0, delta, 2 delta, ...`, linear in
#' between and 0 beyond, `par = delta` (as [jtcal()] tabulates).
#'
#' @param d Distances: numbers from 0.
#' @param type The form: one of the names above.
#' @param par The parameters of the form, of the length given above.
#' @param tab For `type = "table"`: at least two values from 0.
#' @return A numeric vector of `f(d)`, one per distance.
#' @references Brantingham, P. J. and Brantingham, P. L. (1981).
#'   *Environmental Criminology*. Sage.
#'
#'   Levine, N. (2009). Introduction to the special issue on Bayesian
#'   journey-to-crime modelling. *Journal of Investigative Psychology and
#'   Offender Profiling*, 6(3), 167-185. \doi{10.1002/jip.107}
#' @examples
#' jtdecay(c(0, 0.5, 1, 2, 5), type = "trunc", par = c(1, 0.5, 0.8))
#' @export
jtdecay <- function(d, type, par, tab = NULL) {
  .Call(C_clab_r_jtdecay, d, type, par, tab)
}

#' Journey-to-crime surface
#'
#' Ranks the cells of a grid as candidate bases of an offender from the
#' locations of a series of incidents. `mode = "sum"` gives the classic
#' aggregate surface, the sum over incidents of the decay of their distance
#' to the cell (Canter et al. 2000; Levine 2009); `mode = "lik"` gives the
#' sum of the logarithms, the log-likelihood of the cell as base when the
#' decay function is the planar density of an incident about the base and
#' the incidents are independent (O'Leary 2009), so that `prob` is the
#' posterior over cells under a uniform prior. The sum surface is not a
#' likelihood.
#'
#' What this does not license: the surface ranks places to search under
#' the assumptions that the series is one offender's, that the offender
#' travelled from one fixed base inside the grid, and that the decay
#' function, which describes offenders in general, fits this one. It is a
#' search priority, not evidence that any person lives anywhere. Accuracy on
#' solved series in the same place ([jtacc()]) is the only guide to its
#' worth.
#'
#' @inheritParams nnidx
#' @param win The grid's rectangle `c(xmin, xmax, ymin, ymax)`.
#' @param nx,ny The numbers of cells across and up: whole numbers from 1,
#'   with at most 1e9 cells.
#' @inheritParams jtdecay
#' @param metric `"euclid"` (straight-line distance) or `"manhattan"`
#'   (`|dx| + |dy|`, a street grid aligned with the axes).
#' @param mode `"sum"` or `"lik"`.
#' @return A list. `raw` and `prob`, matrices of `nx` rows (across) and
#'   `ny` columns (up), cell `[i, j]` centred at `(xmin + (i - 1/2) dx,
#'   ymin + (j - 1/2) dy)`: the surface and its normalised version (summing
#'   to 1). `peak`, a named vector: `cell` (the first cell of maximum value,
#'   numbered down the columns of the matrix), `x` and `y` (its centre).
#' @references Canter, D., Coffey, T., Huntley, M. and Missen, C. (2000).
#'   Predicting serial killers' home base using a decision support system.
#'   *Journal of Quantitative Criminology*, 16(4), 457-478.
#'   \doi{10.1023/A:1007551316253}
#'
#'   O'Leary, M. (2009). The mathematics of geographic profiling. *Journal
#'   of Investigative Psychology and Offender Profiling*, 6(3), 253-265.
#'   \doi{10.1002/jip.111}
#' @examples
#' s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
#'             ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
#' s$peak
#' @export
jtsurf <- function(x, y = NULL, win, nx, ny, type, par, tab = NULL,
                   metric = "euclid", mode = "sum") {
  .Call(C_clab_r_jtsurf, x, y, win, nx, ny, type, par, tab, metric, mode)
}

#' Empirical distance decay from solved trips
#'
#' A kernel density estimate of trip lengths from known offender trips,
#' to calibrate a decay function (Levine 2009; Silverman 1986). With
#' `refl = TRUE` the kernels are reflected at 0, so no mass falls below 0
#' and the estimate integrates to 1 over `[0, Inf)`. Tabulated at evenly
#' spaced distances it is a `"table"` decay for [jtdecay()] and [jtsurf()].
#'
#' @param d Trip lengths: numbers from 0.
#' @param t The distances to evaluate at (from 0 when `refl = TRUE`).
#' @param h The bandwidth: one positive number.
#' @param w `NULL` for equal weights, or one weight from 0 per trip, with a
#'   positive sum.
#' @param kern The kernel: `"normal"`, `"quartic"`, `"triangular"`,
#'   `"uniform"` or `"laplace"` (double exponential).
#' @param refl `TRUE` to reflect at 0.
#' @return A numeric vector of densities, one per value of `t`.
#' @references Silverman, B. W. (1986). *Density Estimation for Statistics
#'   and Data Analysis*. Chapman and Hall. \doi{10.1201/9781315140919}
#'
#'   Levine, N. (2009). Introduction to the special issue on Bayesian
#'   journey-to-crime modelling. *Journal of Investigative Psychology and
#'   Offender Profiling*, 6(3), 167-185. \doi{10.1002/jip.107}
#' @examples
#' set.seed(1)
#' trips <- rexp(200, 1 / 2)
#' jtcal(trips, t = seq(0, 10, 0.5), h = 0.5)
#' @export
jtcal <- function(d, t, h, w = NULL, kern = "normal", refl = TRUE) {
  .Call(C_clab_r_jtcal, d, t, h, w, kern, refl)
}

#' Accuracy of a journey-to-crime surface on a solved series
#'
#' Given a surface over the grid of [jtsurf()] and the offender's known
#' base: the surface's value in the base's cell, the hit score (the share
#' of the grid searched, in decreasing order of the surface, up to and
#' including the base's cell), in the best and the worst order of ties
#' (Rossmo 2000; Canter et al. 2000), and the error distance from the peak
#' cell's centre to the base.
#'
#' @param prob The surface: a matrix of values from 0, `nx` rows by `ny`
#'   columns, as [jtsurf()] returns.
#' @param win The grid's rectangle, as given to [jtsurf()].
#' @param base The known base, `c(x, y)`, inside `win`.
#' @return A named numeric vector: `p` (the value in the base's cell),
#'   `hitbest` and `hitworst` (shares from 0 to 1), `err` (the error
#'   distance), `peak` and `cell` (the peak's and the base's cell, numbered
#'   down the columns).
#' @references Rossmo, D. K. (2000). *Geographic Profiling*. CRC Press.
#' @examples
#' s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
#'             ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
#' jtacc(s$prob, win = c(0, 6, 0, 6), base = c(3.2, 2.9))
#' @export
jtacc <- function(prob, win, base) {
  .Call(C_clab_r_jtacc, prob, win, base)
}

#' Bayesian journey to crime
#'
#' Combines a journey-to-crime surface over origin zones with an
#' origin-destination matrix of solved offences (Levine 2009; Levine and
#' Block 2011): `P(Jtc)`, the surface; `P(O)`, the origins of all solved
#' offences; `P(O|Jtc)`, the origins of offences in the zones where the
#' series happened; their product, Bayes' numerator with the surface as
#' prior; and the risk, the product over `P(O)`, renormalised. Each block
#' sums to 1. The matrix describes past solved offences and inherits their
#' selection.
#'
#' @param pj The journey-to-crime value of each origin zone: numbers from
#'   0, not all 0.
#' @param od The origin-destination counts: a matrix with one row per origin
#'   zone and one column per destination zone, numbers from 0.
#' @param dest The destination zone (column of `od`) of each incident of
#'   the series: whole numbers from 1; repeats count.
#' @return A data frame with one row per origin zone: `pj`, `po`, `pod`,
#'   `prod` and `risk` (0 where `po` is 0).
#' @references Levine, N. and Block, R. (2011). Bayesian journey-to-crime
#'   estimation: an improvement in geographic profiling methodology. *The
#'   Professional Geographer*, 63(2), 213-229.
#'   \doi{10.1080/00330124.2010.547152}
#' @examples
#' od <- matrix(c(5, 1, 0, 2, 4, 1, 0, 1, 6), 3)
#' bjtc(pj = c(0.5, 0.3, 0.2), od = od, dest = c(1, 2, 2))
#' @export
bjtc <- function(pj, od, dest) {
  .Call(C_clab_r_bjtc, pj, od, dest)
}
