#' Inverse-distance-weighted interpolation, and zonal means
#'
#' Carries the values of monitors to target points as the mean of the
#' nearest monitors weighted by their distance to the power `-p`, and
#' averages the targets of each zone. Exposure papers assign monitor
#' readings to places this way: Rahman et al. (2023, sec. 2.3) take
#' "inverse distance squared weighting of observations from up to four
#' nearby (<50 km) monitoring locations" for census tract centroids (`p =
#' 2`, `k = 4`, `rad = 50`, `lonlat = TRUE`); Wongchanapai et al. (2025)
#' interpolate a surface of annual PM2.5 "using the inverse distance weight
#' method" and take "zonal statistics" per municipality (targets on a grid
#' over the state, `zone` the municipality of each grid point; they used
#' ArcGIS Pro and state neither the power nor the neighbours, so give
#' both); Wesselbaum (2022) matches each crime "to the measures of the
#' nearest weather, pollution station" (`k = 1`).
#'
#' A target at a monitor's location takes that monitor's value (the mean
#' of the monitors there); else the `k` nearest monitors within `rad`
#' (`k = 0`: all within `rad`), ties in distance broken by the order of the
#' monitors, weigh `d^-p` each; `p = 0` gives their plain mean. Monitors
#' with a missing value are left out. A target with no monitor within
#' `rad` is `NA`. With `zone`, the zonal mean is the mean of the zone's
#' targets that have a value: with targets on a regular grid, the mean of
#' the interpolated surface over the zone.
#'
#' @param mx,my The coordinates of the monitors.
#' @param mv The value of each monitor (`NA` allowed).
#' @param tx,ty The coordinates of the targets.
#' @param p The power of the distance, from 0.
#' @param k The number of nearest monitors used, a whole number from 0
#'   (0: all within `rad`).
#' @param rad The radius: a positive number or `Inf`.
#' @param lonlat `TRUE` for longitude and latitude in degrees (distances
#'   in kilometres on the sphere); else planar.
#' @param zone `NULL` or the zone of each target (labels).
#' @return Without `zone`, a data frame, one row per target: `val` and
#'   `used`, the monitors used. With `zone`, one row per zone (sorted):
#'   `zone`, `val` (the zonal mean) and `n`, the targets with a value.
#' @references
#'   Rahman, M. M., Lorenzo, M., Ban-Weiss, G., Hasan, Z., Azzouz, M.,
#'   Eckel, S. P., Conti, D. V., Lurmann, F., Schlaerth, H., Johnston, J.,
#'   Ko, J., Palinkas, L., Hurlburt, M., Silva, S., Gauderman, W. J.,
#'   McConnell, R. and Garcia, E. (2023). Ambient temperature and air
#'   pollution associations with suicide and homicide mortality in
#'   California: a statewide case-crossover study. *Science of the Total
#'   Environment*, 874, 162462. \doi{10.1016/j.scitotenv.2023.162462}
#'
#'   Wongchanapai, P., Green, R., Ku, C. S., Peltier, M. R. and Ranjan, S.
#'   (2025). Air pollution and intimate partner violence in New Jersey:
#'   examining the evidence. *Psychology of Violence*, 15(6), 682-688.
#'   \doi{10.1037/vio0000653}
#'
#'   Wesselbaum, D. (2022). Violent crimes and homicide in New York City:
#'   the role of weather and pollution. *Journal of Forensic and Legal
#'   Medicine*, 91, 102430. \doi{10.1016/j.jflm.2022.102430}
#' @seealso [kerexp()], [bufexp()], [ccross()]
#' @examples
#' m <- data.frame(x = c(0, 10, 0), y = c(0, 0, 10), pm = c(8, 12, 10))
#' idwint(m$x, m$y, m$pm, c(5, 0), c(0, 0))
#' # the nearest monitor only
#' idwint(m$x, m$y, m$pm, c(9, 1), c(0, 1), k = 1)
#' # a grid over two zones, and the zonal means
#' g <- expand.grid(x = 0:10, y = 0:10)
#' idwint(m$x, m$y, m$pm, g$x, g$y, zone = ifelse(g$x < 5, "west", "east"))
#' @export
idwint <- function(mx, my, mv, tx, ty, p = 2, k = 0, rad = Inf,
                   lonlat = FALSE, zone = NULL) {
  f <- "idwint"
  a <- .idwchk(mx, my, mv, tx, ty, p, k, rad, lonlat, f)
  p <- a$p
  k <- a$k
  lonlat <- a$lonlat
  g <- .idwzn(zone, tx, f)
  r <- .Call(C_clab_r_idwint, as.numeric(mx), as.numeric(my), as.numeric(mv),
             as.numeric(tx), as.numeric(ty), p, as.integer(k), as.numeric(rad),
             lonlat, if (is.null(g)) NULL else match(zone, g) - 1L, length(g))
  for (j in seq_along(r)) r[[j]][is.na(r[[j]])] <- NA_real_
  if (is.null(g)) return(data.frame(val = r[[1L]], used = r[[2L]]))
  data.frame(zone = g, val = r[[3L]], n = r[[4L]])
}

# The checks of idwint() in its order: the monitors, the targets, p, k,
# rad and lonlat; list(p, k, lonlat) as checked.
.idwchk <- function(mx, my, mv, tx, ty, p, k, rad, lonlat, f) {
  .cknum(mx, "mx", f, scalar = FALSE)
  .cknum(my, "my", f, scalar = FALSE)
  if (!is.numeric(mv) || length(my) != length(mx) ||
        length(mv) != length(mx)) {
    .ckstop(f, "`mx`, `my` and `mv` must be numeric vectors of one length")
  }
  .cknum(tx, "tx", f, scalar = FALSE)
  .cknum(ty, "ty", f, scalar = FALSE)
  if (length(ty) != length(tx)) .ckstop(f, "`tx` and `ty` must have one length")
  p <- .cknum(p, "p", f, min = 0)
  k <- .cknum(k, "k", f, min = 0, max = 2147483647, whole = TRUE)
  if (!is.numeric(rad) || length(rad) != 1L || is.na(rad) || !(rad > 0)) {
    .ckstop(f, "`rad` must be a positive number or Inf")
  }
  list(p = p, k = k, lonlat = .ckflag(lonlat, "lonlat", f))
}

# The zones of idwint()'s targets, sorted (NULL without zone).
.idwzn <- function(zone, tx, f) {
  if (is.null(zone)) return(NULL)
  if (length(zone) != length(tx) || anyNA(zone)) {
    .ckstop(f, "`zone` must give the zone of each target, with no NA")
  }
  sort(unique(zone))
}

#' Geographically and temporally weighted regression
#'
#' The regression of Kuo and Putra (2021, eqs 4-9), after Huang, Wu and
#' Barry (2010): at each observation `i` the least-squares fit of all
#' observations weighted by a Gaussian kernel of their space-time distance,
#' so that the coefficients vary over places and periods. They relate
#' robbery, domestic violence and fraud in the local government areas of
#' New South Wales to the air quality index month by month.
#'
#' With the spatial distance `ds` (planar, or kilometres on the sphere with
#' `lonlat`) and the time distance `dt`, the space-time distance is
#' `d^2 = ds^2 + tau dt^2` (eq 7, the spatial scale set to 1) and the
#' weight `w = exp(-d^2 / h^2)` (eq 8); `beta_i = (X'W_i X)^-1 X'W_i y`
#' (eq 5). Without `t` (or with `tau = 0`) it is geographically weighted
#' regression. The cross-validation score is `sum (y_i - y_(-i))^2` (eq
#' 9), the leave-one-out prediction from the fit at `i` without
#' observation `i` (computed as `e_i / (1 - S_ii)`, `S` the hat matrix).
#' Kuo and Putra choose `tau` and the bandwidth `h` by cross-validation or
#' by "the corrected Akaike information criterion" without giving it; the
#' `aicc` here is the one of the GWR literature, `n log(RSS / n) + n log(2
#' pi) + n (n + tr S) / (n - 2 - tr S)`. Give several values of `h` and
#' `tau` and the pair with the smallest `crit` is fitted; all are in
#' `grid`. The local standard errors are `sigma^2` times the diagonal of
#' `(X'W_iX)^-1 X'W_i^2 X (X'W_iX)^-1`, `sigma^2 = RSS / (n - 2 tr S +
#' tr S'S)`.
#'
#' @param y The outcome: finite numbers.
#' @param x The regressors: a numeric matrix or data frame, one row per
#'   value of `y` (columns named), or `NULL` for the constant alone.
#' @param u,v The coordinates of each observation.
#' @param t `NULL` or the time of each observation (numbers).
#' @param h The bandwidths tried: positive numbers.
#' @param tau The ratios of the time to the space scale tried: numbers
#'   from 0.
#' @param crit `"cv"` or `"aicc"`: how `h` and `tau` are chosen.
#' @param lonlat `TRUE` for longitude and latitude in degrees.
#' @param cons `TRUE` for a local constant.
#' @return An object of class `gtwreg`, a list: `coef` and `se`, matrices
#'   of the local coefficients and their standard errors (one row per
#'   observation); `fit`, `resid`, `hat` and `loo` (the leave-one-out
#'   residuals); `stat`, a named vector (`h`, `tau`, `rss`, `trs`, `trss`,
#'   `sigma2`, `cv`, `aicc`, `edf`, `r2`); and `grid`, the criteria of every
#'   pair tried.
#' @references
#'   Kuo, P.-F. and Putra, I. G. B. (2021). Analyzing the relationship
#'   between air pollution and various types of crime. *PLoS ONE*, 16(8),
#'   e0255653. \doi{10.1371/journal.pone.0255653}
#'
#'   Huang, B., Wu, B. and Barry, M. (2010). Geographically and temporally
#'   weighted regression for modeling spatio-temporal variation in house
#'   prices. *International Journal of Geographical Information Science*,
#'   24(3), 383-401. \doi{10.1080/13658810802672469}
#' @seealso [sarml()], [gmoran()]
#' @examples
#' set.seed(4)
#' d <- expand.grid(u = 1:6, v = 1:6, t = 1:4)
#' d$aqi <- rnorm(nrow(d))
#' d$y <- 1 + (0.2 + 0.1 * d$u) * d$aqi + rnorm(nrow(d), 0, 0.3)
#' r <- gtwreg(d$y, d["aqi"], d$u, d$v, d$t, h = c(1.5, 3, 6),
#'             tau = c(0, 1))
#' r
#' r$grid
#' @export
gtwreg <- function(y, x, u, v, t = NULL, h, tau = 0, crit = "cv",
                   lonlat = FALSE, cons = TRUE) {
  f <- "gtwreg"
  a <- .gtchk(y, x, u, v, t, cons, f)
  h <- .cknum(h, "h", f, scalar = FALSE, positive = TRUE)
  tau <- .cknum(tau, "tau", f, scalar = FALSE, min = 0)
  crit <- .ckarg(crit, c("cv", "aicc"), "crit", f)
  lonlat <- .ckflag(lonlat, "lonlat", f)
  g <- expand.grid(h = h, tau = tau)
  run <- function(hh, tt) {
    .Call(C_clab_r_gtwreg, a$x, a$y, a$u, a$v, a$t, hh, tt, lonlat)
  }
  g$cv <- NA_real_
  g$aicc <- NA_real_
  for (i in seq_len(nrow(g))) {
    s <- run(g$h[i], g$tau[i])$stat
    g$cv[i] <- s[["cv"]]
    g$aicc[i] <- s[["aicc"]]
  }
  b <- which.min(g[[crit]])
  if (!length(b)) .ckstop(f, "no bandwidth gives a finite `%s`", crit)
  r <- run(g$h[b], g$tau[b])
  .gtout(r, a, g, g$h[b], g$tau[b])
}

#' @export
print.gtwreg <- function(x, ...) {
  s <- x[["stat"]]
  cat("gtwreg: geographically and temporally weighted regression, h = ",
      format(s[["h"]], digits = 4), ", tau = ", format(s[["tau"]], digits = 4),
      "\nCV ", format(s[["cv"]], digits = 5), ", AICc ",
      format(s[["aicc"]], digits = 5), ", R2 ", format(s[["r2"]], digits = 3),
      ", effective parameters ", format(s[["trs"]], digits = 4),
      "\nlocal coefficients (quartiles):\n", sep = "")
  print(t(apply(x[["coef"]], 2, stats::quantile, c(0, 0.25, 0.5, 0.75, 1))))
  invisible(x)
}

#' Association rules between two categorical variables
#'
#' The rules the Apriori algorithm (Agrawal et al. 1996) finds when each
#' observation holds one level of each of two variables: `{} => {b}` and
#' `{a} => {b}` with `a` a level of one variable and `b` of the other, each
#' with its support, confidence, coverage and lift. Wongchanapai et al.
#' (2025) cut the rate of intimate partner violence and the annual PM2.5 of
#' New Jersey's municipalities into five groups each with [cut()] and kept
#' rules of "at least 10% support ... and at least 20% confidence"
#' (`minsup = 0.1`, `minconf = 0.2`, the defaults here) in R's arules
#' package; with two variables these are all its rules, the empty
#' left-hand side included as arules includes it by default (`empty`).
#'
#' Support is the share of observations holding both sides; confidence,
#' the support over that of the left-hand side (its coverage); lift, the
#' confidence over the support of the right-hand side, 1 for a side that
#' tells nothing of the other. Two levels of one variable never occur
#' together, so no rule joins them, and a pair of levels never observed
#' together is no rule.
#'
#' @param x,y The two variables: factors (their levels, in order) or
#'   vectors (their sorted values), one value per observation, no `NA`.
#' @param minsup,minconf The least support and confidence, from 0 to 1.
#' @param empty `TRUE` to keep the rules with an empty left-hand side.
#' @param nm The names of the two variables in the items.
#' @return A data frame, one row per rule: `lhs`, `rhs` (items `"name=level"`,
#'   `""` for the empty set), `support`, `confidence`, `coverage`, `lift`
#'   and `count`.
#' @references
#'   Agrawal, R., Mannila, H., Srikant, R., Toivonen, H. and Verkamo, A. I.
#'   (1996). Fast discovery of association rules. In *Advances in Knowledge
#'   Discovery and Data Mining*, 307-328. AAAI Press.
#'
#'   Wongchanapai, P., Green, R., Ku, C. S., Peltier, M. R. and Ranjan, S.
#'   (2025). Air pollution and intimate partner violence in New Jersey:
#'   examining the evidence. *Psychology of Violence*, 15(6), 682-688.
#'   \doi{10.1037/vio0000653}
#' @seealso [ctabx()]
#' @examples
#' set.seed(5)
#' pm <- cut(rnorm(200), 5, labels = c("vl", "l", "m", "mh", "vh"))
#' ipv <- cut(rnorm(200), 5, labels = c("vl", "l", "m", "mh", "vh"))
#' asrule(pm, ipv, nm = c("pm25", "ipv"))
#' @export
asrule <- function(x, y, minsup = 0.1, minconf = 0.2, empty = TRUE,
                   nm = c("x", "y")) {
  f <- "asrule"
  a <- .archk(x, y, minsup, minconf, empty, nm, f)
  lx <- .arlev(x)
  ly <- .arlev(y)
  r <- .Call(C_clab_r_asrule, match(as.character(x), as.character(lx)) - 1L,
             length(lx), match(as.character(y), as.character(ly)) - 1L,
             length(ly), a$minsup, a$minconf, a$empty)
  it <- c(paste0(nm[1L], "=", lx), paste0(nm[2L], "=", ly))
  # the empty left-hand side is -1: index 0 would drop out of the labels
  lh <- it[pmax(r$lhs, 0L) + 1L]
  lh[r$lhs < 0L] <- ""
  data.frame(lhs = lh, rhs = it[r$rhs + 1L], support = r$out[, 1L],
             confidence = r$out[, 2L], coverage = r$out[, 3L],
             lift = r$out[, 4L], count = r$out[, 5L])
}

# The checks of asrule() in its order; list(minsup, minconf, empty).
.archk <- function(x, y, minsup, minconf, empty, nm, f) {
  if (length(x) != length(y) || !length(x) || anyNA(x) || anyNA(y)) {
    .ckstop(f, "`x` and `y` must have one length and no NA")
  }
  minsup <- .cknum(minsup, "minsup", f, min = 0, max = 1)
  minconf <- .cknum(minconf, "minconf", f, min = 0, max = 1)
  empty <- .ckflag(empty, "empty", f)
  if (!is.character(nm) || length(nm) != 2L || anyNA(nm)) {
    .ckstop(f, "`nm` must be two names")
  }
  list(minsup = minsup, minconf = minconf, empty = empty)
}

# The items of one column of asrule(): a factor's levels, else the sorted
# values.
.arlev <- function(x) if (is.factor(x)) levels(x) else sort(unique(x))

# The checks of gtwreg() in its order: the outcome, the regressors (with
# the constant first), the coordinates and the times, as doubles.
.gtchk <- function(y, x, u, v, t, cons, f) {
  .cknum(y, "y", f, scalar = FALSE)
  n <- length(y)
  cons <- .ckflag(cons, "cons", f)
  xm <- if (is.null(x)) matrix(0, n, 0) else as.matrix(x)
  if (!is.numeric(xm) || nrow(xm) != n || anyNA(xm) || !all(is.finite(xm))) {
    .ckstop(f, "`x` must be finite numbers, one row per value of `y`")
  }
  if (is.null(colnames(xm))) colnames(xm) <- paste0("x", seq_len(ncol(xm)))
  if (cons) xm <- cbind("(const)" = 1, xm)
  if (!ncol(xm)) .ckstop(f, "the model needs `x` or `cons = TRUE`")
  .cknum(u, "u", f, scalar = FALSE)
  .cknum(v, "v", f, scalar = FALSE)
  if (length(u) != n || length(v) != n) {
    .ckstop(f, "`u` and `v` must have one value per value of `y`")
  }
  if (!is.null(t)) {
    .cknum(t, "t", f, scalar = FALSE)
    if (length(t) != n) .ckstop(f, "`t` must have one value per value of `y`")
  }
  if (n <= ncol(xm) + 1L) {
    .ckstop(f, "the model needs more observations than coefficients plus one")
  }
  list(y = as.numeric(y), x = xm, u = as.numeric(u), v = as.numeric(v),
       t = if (is.null(t)) NULL else as.numeric(t))
}

# The result of gtwreg() from the kernel's r at the chosen h and tau.
.gtout <- function(r, a, g, h, tau) {
  cn <- colnames(a$x)
  dimnames(r$coef) <- dimnames(r$se) <- list(NULL, cn)
  r$loo[is.na(r$loo)] <- NA_real_
  s <- r$stat
  s <- c(h = h, tau = tau, s,
         r2 = 1 - s[["rss"]] / sum((a$y - mean(a$y))^2))
  structure(list(coef = r$coef, se = r$se, fit = r$fit, resid = a$y - r$fit,
                 hat = r$hat, loo = r$loo, stat = s, grid = g),
            class = "gtwreg")
}
