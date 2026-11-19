#' Empirical semivariogram
#'
#' The classical estimator of Matheron (1963): half the mean squared
#' difference of the values over the pairs of points whose distance falls
#' in each of `nbins` equal bins of `(0, cutoff]` (the first bin closed at
#' 0). Bins without pairs are left out. The semivariogram describes how
#' dissimilar values become with distance, the first step of kriging:
#' rates of a crime, a disease or a survey outcome measured at places.
#'
#' @inheritParams cgstat
#' @param v The values: finite numbers, one per point.
#' @param nbins The number of distance bins: one whole number from 1 to
#'   67108864 (2^26).
#' @param cutoff The largest distance used: `NULL` for a third of the
#'   largest distance between two points (the convention of gstat), or one
#'   positive number.
#' @return A data frame with one row per bin holding pairs: `dist` (the
#'   centre of the bin), `gamma` (the semivariance) and `np` (the number of
#'   pairs).
#' @references Matheron, G. (1963). Principles of geostatistics.
#'   *Economic Geology*, 58(8), 1246-1266. \doi{10.2113/gsecongeo.58.8.1246}
#' @examples
#' # simulated values with a trend from west to east
#' set.seed(1)
#' p <- cbind(x = runif(60, 0, 10), y = runif(60, 0, 10))
#' v <- p[, 1] * 0.5 + rnorm(60, 0, 0.3)
#' vgemp(p, v = v)
#' @family geostatistics
#' @export
vgemp <- function(x, y = NULL, v, nbins = 15, cutoff = NULL) {
  .Call(C_clab_r_vgemp, x, y, v, nbins, cutoff)
}

#' Fit a variogram model by Gaussian maximum likelihood
#'
#' Fits the nugget, partial sill and range of a spherical, exponential or
#' Gaussian covariance model by maximising the Gaussian log-likelihood of
#' the values about a constant mean (Mardia and Marshall 1984). The search
#' is the Nelder-Mead simplex of R's `optim()` (at most 800 evaluations)
#' on the logarithms of the three parameters, from 0.1 and 0.9 of the
#' variance of the values and a quarter of the largest distance. The
#' covariance matrix gets `1e-8` of the variance on its diagonal so that
#' it stays positive definite.
#'
#' @inheritParams vgemp
#' @param model `"exponential"`, `"spherical"` or `"gaussian"`.
#' @return A list of class `vgfit`: `model`, `nugget`, `psill`, `range`,
#'   `loglik` (the log-likelihood less `n / 2 log(2 pi)`), `conv` (`TRUE`
#'   when the simplex converged) and `nfn` (the evaluations).
#' @references Mardia, K. V. and Marshall, R. J. (1984). Maximum likelihood
#'   estimation of models for residual covariance in spatial regression.
#'   *Biometrika*, 71(1), 135-146. \doi{10.1093/biomet/71.1.135}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(40, 0, 10), y = runif(40, 0, 10))
#' v <- p[, 1] * 0.5 + rnorm(40, 0, 0.3)
#' vgfit(p, v = v)
#' @family geostatistics
#' @export
vgfit <- function(x, y = NULL, v, model = "exponential") {
  .Call(C_clab_r_vgfit, x, y, v, model)
}

#' @export
print.vgfit <- function(x, ...) {
  cat(
    "vgfit: ",
    x[["model"]],
    " model by maximum likelihood",
    if (!x[["conv"]]) " (not converged)",
    "\n",
    sep = ""
  )
  print(unlist(x[c("nugget", "psill", "range", "loglik")]))
  invisible(x)
}

#' Ordinary kriging
#'
#' Predicts the values at new locations `at` by ordinary kriging: the
#' best linear unbiased predictor under a variogram model with an unknown
#' constant mean (Cressie 1993, sec. 3.2), with its kriging variance. The
#' variogram is `vgm`, a list such as [vgfit()] returns; without it, an
#' exponential model is fitted by weighted least squares to the 15-bin
#' [vgemp()] semivariogram with the weights `np / h^2` of Cressie (1985),
#' which is fast but rougher than [vgfit()].
#'
#' @inheritParams vgemp
#' @param at The prediction locations: a matrix or data frame of two
#'   coordinate columns, as `x`.
#' @param vgm `NULL`, or a list with `model` (`"exponential"`,
#'   `"spherical"` or `"gaussian"`), `nugget`, `psill` (from 0) and `range`
#'   (positive).
#' @return A data frame with one row per prediction location: `pred` and
#'   `var` (the kriging variance, at least 0).
#' @references Cressie, N. (1985). Fitting variogram models by weighted
#'   least squares. *Mathematical Geology*, 17(5), 563-586.
#'   \doi{10.1007/BF01032109}
#'
#'   Cressie, N. (1993). *Statistics for Spatial Data*, revised edition.
#'   Wiley, New York. \doi{10.1002/9781119115151}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(60), y = runif(60))
#' v <- sin(3 * p[, 1]) + rnorm(60, sd = 0.1)
#' okrge(p, v = v, at = cbind(0.5, 0.5))
#' okrge(p, v = v, at = cbind(0.5, 0.5),
#'       vgm = list(model = "exponential", nugget = 0.01, psill = 0.3,
#'                  range = 0.4))
#' @family geostatistics
#' @export
okrge <- function(x, y = NULL, v, at, vgm = NULL) {
  .Call(C_clab_r_okrge, x, y, v, at, vgm)
}

#' Ordinary block kriging
#'
#' Predicts the mean value over each block (an area such as a district)
#' by ordinary kriging with point-to-block covariances averaged over the
#' block, `C(s_i, B) = (1 / |B|) int_B` `C(s_i, u) du`, under the exponential
#' covariance `(sill - nugget) exp(-h / range)`, plus `nugget` at distance
#' 0 (Schabenberger and Gotway 2005, ch. 4). A block is a box, given as a
#' 2 x 2 matrix of its lower and upper corners (rows) and integrated on a
#' regular grid of `round(sqrt(nq))^2` points, or a matrix of the points
#' that discretise it.
#'
#' @inheritParams vgemp
#' @param blocks A list of blocks: each a 2 x 2 matrix (rows: lower and
#'   upper corner, columns: x and y), or a matrix of two columns of points.
#' @param nq The size of the grid of a box: one whole number from 1. A box
#'   is averaged over `round(sqrt(nq))^2` points (`nq` rounded to a square:
#'   10 gives 9). The points of all the blocks must be at most 67108864
#'   (2^26).
#' @param nugget,sill,range The covariance: `nugget` from 0, `sill` at least
#'   `nugget`, `range` positive.
#' @return A data frame with one row per block: `est` (the kriged mean) and
#'   `se` (its kriging standard error).
#' @references Schabenberger, O. and Gotway, C. A. (2005). *Statistical
#'   Methods for Spatial Data Analysis*. Chapman and Hall/CRC, Boca Raton.
#'   \doi{10.1201/9781315275086}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(30), y = runif(30))
#' v <- p[, 1] + rnorm(30, sd = 0.1)
#' bkrig(p, v = v, blocks = list(rbind(c(0, 0), c(0.5, 0.5)),
#'                               rbind(c(0.5, 0.5), c(1, 1))),
#'       range = 0.3)
#' @family geostatistics
#' @export
bkrig <- function(
  x,
  y = NULL,
  v,
  blocks,
  nq = 25,
  nugget = 0,
  sill = 1,
  range = 1
) {
  .Call(C_clab_r_bkrig, x, y, v, blocks, nq, nugget, sill, range)
}

#' Leave-one-out cross-validation of ordinary kriging
#'
#' Predicts each value from all the others by ordinary kriging with the
#' exponential covariance of [bkrig()], and summarises the residuals by
#' the mean squared prediction error `MSPE = (1 / n) sum (Z(s_i) -`
#' `Z_-i(s_i))^2`, its root and the mean absolute error (Schabenberger and
#' Gotway 2005, ch. 4). Comparing it across covariance parameters chooses
#' among them.
#'
#' @inheritParams bkrig
#' @return A list of class `krigcv`: `mspe`, `rmspe`, `mae` and `resid` (the
#'   leave-one-out residuals).
#' @inherit bkrig references
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(30), y = runif(30))
#' krigcv(p, v = rnorm(30), range = 0.2)
#' @family geostatistics
#' @export
krigcv <- function(x, y = NULL, v, nugget = 0, sill = 1, range = 1) {
  .Call(C_clab_r_krigcv, x, y, v, nugget, sill, range)
}

#' @export
print.krigcv <- function(x, ...) {
  cat(
    "krigcv: leave-one-out kriging of ",
    length(x[["resid"]]),
    " values\n",
    sep = ""
  )
  print(unlist(x[c("mspe", "rmspe", "mae")]))
  invisible(x)
}

#' Polynomial trend surface
#'
#' Fits the mean of the values as a polynomial in the coordinates by
#' ordinary least squares, `mu(s) = sum_k beta_k f_k(s)`, the `f_k` the
#' monomials up to degree `order`: `1`; `x`, `y`; `x2`, `y2`, `xy`; `x3`,
#' `y3`, `x2y`, `xy2` (Schabenberger and Gotway 2005, ch. 2). The large-
#' scale trend is often removed this way before a variogram is fitted to
#' the residuals.
#'
#' @inheritParams vgemp
#' @param order The degree: 0, 1, 2 or 3.
#' @return A list of class `trsurf`: `coef` and `se` (named by monomial),
#'   `r2`, `order` and `n`. [summary()] gives the table of coefficients.
#' @inherit bkrig references
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(50), y = runif(50))
#' v <- 1 + 2 * p[, 1] - p[, 2]^2 + rnorm(50, sd = 0.1)
#' f <- trsurf(p, v = v)
#' f
#' summary(f)
#' @family geostatistics
#' @export
trsurf <- function(x, y = NULL, v, order = 2) {
  .Call(C_clab_r_trsurf, x, y, v, order)
}

#' @export
print.trsurf <- function(x, ...) {
  cat(
    "trsurf: trend surface of order ",
    x[["order"]],
    " on ",
    x[["n"]],
    " points, R-squared ",
    format(signif(x[["r2"]], 4)),
    "\n",
    sep = ""
  )
  print(x[["coef"]])
  invisible(x)
}

#' @export
summary.trsurf <- function(object, ...) .ctab(object[["coef"]], object[["se"]])

#' Spatial linear model by profile maximum likelihood
#'
#' The Gaussian spatial linear model `v = Z beta + delta + e`, `delta` a
#' Gaussian process with covariance `sigma2 exp(-h / phi)`, without a
#' nugget (Schabenberger and Gotway 2005, ch. 5). The range `phi` is found
#' by maximising the profile likelihood with Brent's method on `log(phi)`
#' over `[log(max(h / 100, 0.001)),` `log(max(3 h, 1))]`, `h` the largest
#' distance, as R's `optimize()` runs it; `beta` follows by generalised
#' least squares. The standard errors use `sigma2` estimated with `n - p`
#' degrees of freedom.
#'
#' @inheritParams vgemp
#' @param v The outcome: finite numbers, one per point.
#' @param z The design: `NULL` for an intercept alone, or a numeric matrix
#'   (or vector, or data frame) with one row per point; an intercept is
#'   not added.
#' @return A list of class `sglmx`: `coef`, `se`, `sigma2`, `phi`, `tau2`
#'   (the nugget, 0) and `n`. [summary()] gives the table of coefficients.
#' @inherit bkrig references
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(40), y = runif(40))
#' z <- cbind(1, p[, 1])
#' sglmx(1 + 2 * p[, 1] + rnorm(40, sd = 0.2), z, p)
#' @family geostatistics
#' @export
sglmx <- function(v, z = NULL, x, y = NULL) {
  if (is.data.frame(z)) {
    z <- as.matrix(z)
  }
  .Call(C_clab_r_sglmx, v, z, x, y)
}

#' @export
print.sglmx <- function(x, ...) {
  cat(
    "sglmx: spatial linear model, exponential covariance, ",
    x[["n"]],
    " points\n",
    sep = ""
  )
  print(c(x[["coef"]], sigma2 = x[["sigma2"]], phi = x[["phi"]]))
  invisible(x)
}

#' @export
summary.sglmx <- function(object, ...) .ctab(object[["coef"]], object[["se"]])

#' Spatial linear mixed model by REML
#'
#' The model `v = Z beta + delta + e`, `delta ~ N(0, sigma2 R_phi)` with the
#' exponential correlation `R_phi = exp(-h / phi)` and a nugget `e ~ N(0,
#' tau2 I)`, `tau2 = nu sigma2`, by restricted maximum likelihood
#' (Patterson and Thompson 1971; Schabenberger and Gotway 2005, ch. 5). The
#' search is R's Nelder-Mead (`optim()`, 400 evaluations, relative
#' tolerance `1e-6`) on `(log phi, log nu)` from `(log(h / 3), log(0.1))`,
#' `h` the largest distance (at least 1).
#'
#' @inheritParams sglmx
#' @return A list of class `smixm`: `coef`, `se`, `sigma2`, `tau2`, `phi`,
#'   `conv` and `n`. [summary()] gives the table of coefficients.
#' @references Patterson, H. D. and Thompson, R. (1971). Recovery of
#'   inter-block information when block sizes are unequal. *Biometrika*,
#'   58(3), 545-554. \doi{10.1093/biomet/58.3.545}
#'
#'   Schabenberger, O. and Gotway, C. A. (2005). *Statistical Methods for
#'   Spatial Data Analysis*. Chapman and Hall/CRC, Boca Raton.
#'   \doi{10.1201/9781315275086}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(40), y = runif(40))
#' smixm(1 + 2 * p[, 1] + rnorm(40, sd = 0.3), cbind(1, p[, 1]), p)
#' @family geostatistics
#' @export
smixm <- function(v, z = NULL, x, y = NULL) {
  if (is.data.frame(z)) {
    z <- as.matrix(z)
  }
  .Call(C_clab_r_smixm, v, z, x, y)
}

#' @export
print.smixm <- function(x, ...) {
  cat(
    "smixm: spatial linear mixed model by REML, ",
    x[["n"]],
    " points",
    if (!x[["conv"]]) " (not converged)",
    "\n",
    sep = ""
  )
  print(c(
    x[["coef"]],
    sigma2 = x[["sigma2"]],
    tau2 = x[["tau2"]],
    phi = x[["phi"]]
  ))
  invisible(x)
}

#' @export
summary.smixm <- function(object, ...) .ctab(object[["coef"]], object[["se"]])

#' A spatial data object
#'
#' Holds points with their attributes and their coordinate reference
#' system: the coordinates come from two columns of a data frame, the
#' other columns are kept as attributes. The system is an EPSG code, never
#' a PROJ string: 4326 (WGS 84 longitude and latitude, in degrees), 3857
#' (Web Mercator, in metres) or another code, taken to be in metres.
#' [spxy()] returns the coordinates, [crsof()] the system, [lfspec()] the
#' specification of a web map and [spsamp()] a sample of the points.
#'
#' @param data A data frame holding the coordinate columns and any
#'   attributes.
#' @param coords The names of the two coordinate columns: x (or longitude)
#'   first.
#' @param crs The EPSG code: one whole number.
#' @return A list of class `spobj`: `coords` (a two-column matrix), `data`
#'   (the other columns), `cnames`, `crs`, `units` and `n`.
#' @examples
#' d <- data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7), v = c(1, 2))
#' s <- spobj(d, c("lon", "lat"))
#' s
#' format(s)
#' summary(s)
#' @family geostatistics
#' @export
spobj <- function(data, coords = c("lon", "lat"), crs = 4326) {
  .spargs(data, coords, crs)
  xy <- as.matrix(data[coords])
  if (!is.numeric(xy) || !all(is.finite(xy))) {
    stop("spobj: the coordinates must be finite numbers.", call. = FALSE)
  }
  storage.mode(xy) <- "double"
  rownames(xy) <- NULL
  structure(
    list(
      coords = xy,
      data = data[setdiff(names(data), coords)],
      cnames = coords,
      crs = crs,
      units = if (crs == 4326) "degrees" else "metres",
      n = nrow(xy)
    ),
    class = "spobj"
  )
}

#' @export
format.spobj <- function(x, ...) {
  sprintf(
    "spobj: %d points, EPSG:%s (%s)",
    x[["n"]],
    format(x[["crs"]]),
    x[["units"]]
  )
}

#' @export
print.spobj <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
summary.spobj <- function(object, ...) {
  xy <- object[["coords"]]
  data.frame(
    coord = object[["cnames"]],
    min = apply(xy, 2L, min),
    max = apply(xy, 2L, max),
    row.names = NULL
  )
}

#' The coordinates of a spatial data object
#'
#' @param x A [spobj()] object.
#' @return The two-column matrix of coordinates, named by their columns.
#' @examples
#' spxy(spobj(data.frame(lon = 1:3, lat = 1:3)))
#' @family geostatistics
#' @export
spxy <- function(x) {
  .isobj(x, "spxy")
  x[["coords"]]
}

#' The coordinate reference system of a spatial data object
#'
#' Gives the system as its EPSG code and a label in the style of WKT2,
#' never as a PROJ string.
#'
#' @inheritParams spxy
#' @return A list: `epsg` and `wkt`.
#' @examples
#' crsof(spobj(data.frame(lon = 1, lat = 1)))
#' @family geostatistics
#' @export
crsof <- function(x) {
  .isobj(x, "crsof")
  if (!is.numeric(x[["crs"]]) || length(x[["crs"]]) != 1L) {
    .ckstop("crsof", "`x$crs` must be one EPSG code")
  }
  wkt <- switch(
    as.character(x[["crs"]]),
    "4326" = "GEOGCRS[\"WGS 84\", EPSG:4326]",
    "3857" = "PROJCRS[\"WGS 84 / Pseudo-Mercator\", EPSG:3857]",
    sprintf("EPSG:%s", x[["crs"]])
  )
  list(epsg = x[["crs"]], wkt = wkt)
}

#' The specification of a web map
#'
#' Describes an interactive map of the points as data, for a mapping
#' library to draw: the tile provider, the centre (the mean longitude and
#' latitude), a zoom level and the points in longitude and latitude. Web
#' Mercator coordinates (EPSG:3857) are taken back to longitude and
#' latitude by the inverse of the spherical Mercator projection.
#'
#' @inheritParams spxy
#' @return A list: `provider`, `interactive`, `center` (`lon`, `lat`),
#'   `zoom` and `points` (a matrix of `lon` and `lat`).
#' @examples
#' lfspec(spobj(data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7))))
#' @family geostatistics
#' @export
lfspec <- function(x) {
  .Call(C_clab_r_lfspec, x)
}

#' Sample points, by local density or uniformly
#'
#' Draws `size` of the points without replacement. With `dens = TRUE` each
#' draw picks a point with probability proportional to its local density
#' `sum_j exp(-(d_ij / b)^2)`, `b` the median distance between points
#' (great-circle distance for EPSG:4326), so that dense areas are sampled
#' more; otherwise every point is equally likely.
#'
#' @inheritParams spxy
#' @param size The number of points: one whole number from 0 to the number
#'   of points.
#' @param dens `TRUE` to sample by local density.
#' @param seed `NULL` (a seed drawn from R's generator, so that
#'   `set.seed()` fixes the draws) or one whole number from 0.
#' @return An integer vector of row numbers, with the attribute `prob`, the
#'   probability of each point at the first draw.
#' @examples
#' s <- spobj(expand.grid(lon = 1:5, lat = 1:5))
#' spsamp(s, 5, seed = 1)
#' @family geostatistics
#' @export
spsamp <- function(x, size, dens = TRUE, seed = NULL) {
  .Call(C_clab_r_spsamp, x, size, dens, seed)
}

# the coefficient table of a fit
.ctab <- function(b, s) {
  data.frame(
    term = names(b),
    est = unname(b),
    se = unname(s),
    z = unname(b / s)
  )
}

.spargs <- function(data, coords, crs) {
  if (!is.data.frame(data)) {
    stop("spobj: `data` must be a data frame.", call. = FALSE)
  }
  if (
    !is.character(coords) ||
      length(coords) != 2L ||
      anyNA(coords) ||
      !all(coords %in% names(data))
  ) {
    stop("spobj: `coords` must name two columns of `data`.", call. = FALSE)
  }
  if (
    !is.numeric(crs) ||
      length(crs) != 1L ||
      is.na(crs) ||
      crs != round(crs) ||
      crs < 1
  ) {
    stop("spobj: `crs` must be one EPSG code, a whole number.", call. = FALSE)
  }
}

.isobj <- function(x, f) {
  if (!inherits(x, "spobj")) {
    stop(f, ": `x` must be the object spobj() returns.", call. = FALSE)
  }
}
