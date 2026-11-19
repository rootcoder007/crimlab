#' Thin-plate spline surface
#'
#' `tpsurf()` fits a thin-plate spline to values measured at points in the
#' plane (rates at area centroids, say) and evaluates it on a grid over
#' their box: a smooth surface that passes through the values when `lam`
#' is 0 and is smoothed towards a plane as `lam` grows. It complements
#' [kdsurf()], which smooths the points themselves rather than values at
#' them.
#'
#' @details The surface is `f(s) = sum_i a_i phi(|s -`
#'   `s_i|) + b_0 + b_1 x + b_2 y` with the radial basis `phi(r) = r^2 log r`
#'   (Duchon 1977); the
#'   coefficients solve the linear system `(K + lam I) a + T b = z`, `T'a =
#'   0`, `K` the basis between the points and `T` the rows `(1, x_i, y_i)`
#'   (Wahba 1990, ch. 2). The compiled kernel solves it by LU decomposition
#'   and evaluates the grid. It refuses points that all lie on one line
#'   (the plane is not determined) and, with `lam = 0`, points that repeat
#'   (the system is singular). `fitted` is the surface at the points, so
#'   `resid` is 0 for `lam = 0` up to rounding.
#'
#' @param x,y The coordinates of the points. For `print()` and `plot()`,
#'   `x` is a `tpsurf` result.
#' @param z The value at each point.
#' @param lam The smoothing: one number from 0 (0 interpolates).
#' @param nx,ny The size of the grid: whole numbers from 2 to 1000.
#' @param ... Not used by `print()`; for `plot()`, passed to
#'   [graphics::image()].
#' @inheritParams chomap
#' @param col Not used: the surface is shaded by `pal`.
#' @return An object of class `tpsurf`, a list: `x` and `y` (the grid),
#'   `z` (the surface on it, `nx` by `ny`, as [graphics::image()] takes
#'   it), `a` (the basis coefficients), `beta` (`const`, `x`, `y`),
#'   `fitted`, `resid` and `lam`. `plot()` draws the surface with its
#'   contours and returns `x`, `y` and `z` invisibly.
#' @references Duchon, J. (1977). Splines minimizing rotation-invariant
#'   semi-norms in Sobolev spaces. In W. Schempp and K. Zeller (eds.),
#'   *Constructive Theory of Functions of Several Variables*, Lecture Notes
#'   in Mathematics 571, 85-100. Springer. \doi{10.1007/BFb0086566}
#'
#'   Wahba, G. (1990). *Spline Models for Observational Data*. SIAM.
#'   \doi{10.1137/1.9781611970128}
#' @examples
#' set.seed(2)
#' x <- runif(40)
#' y <- runif(40)
#' s <- tpsurf(x, y, sin(4 * x) + y^2 + rnorm(40, sd = 0.05), lam = 0.001,
#'             nx = 30)
#' s
#' plot(s)
#' @export
tpsurf <- function(x, y, z, lam = 0, nx = 50L, ny = nx) {
  f <- "tpsurf"
  .tpchk(x, y, z, f)
  if (!.pnum1(lam, 0)) .ckstop(f, "`lam` must be one number from 0")
  if (!.pnum1(nx, 2, 1000, TRUE) || !.pnum1(ny, 2, 1000, TRUE)) {
    .ckstop(f, "`nx` and `ny` must be whole numbers from 2 to 1000")
  }
  if (length(x) > 4000L) .ckstop(f, "at most 4000 points")
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  gx <- seq(min(x), max(x), length.out = nx)
  gy <- seq(min(y), max(y), length.out = ny)
  r <- .Call(C_clab_r_tpsurf, x, y, z, as.numeric(lam), gx, gy)
  structure(
    list(x = gx, y = gy, z = r$z, a = r$a,
         beta = c(const = r$b[1L], x = r$b[2L], y = r$b[3L]),
         fitted = r$fit, resid = z - r$fit, lam = as.numeric(lam)),
    class = "tpsurf"
  )
}

#' @rdname tpsurf
#' @export
print.tpsurf <- function(x, ...) {
  cat("tpsurf: thin-plate spline of ", length(x[["a"]]), " points, lambda ",
      format(x[["lam"]]), ", residual sum of squares ",
      format(signif(sum(x[["resid"]]^2), 4)), "; grid ", length(x[["x"]]),
      " x ", length(x[["y"]]), "\n", sep = "")
  invisible(x)
}

#' @rdname tpsurf
#' @export
plot.tpsurf <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "tpsurf", c("x", "y", "z"))
  .psurf(x$x, x$y, x$z, .pmain(main, "Thin-plate spline surface"),
         pal = pal, ...)
}

# The points: numbers of the same length, finite, at least three.
.tpchk <- function(x, y, z, f) {
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(z) ||
        length(y) != length(x) || length(z) != length(x)) {
    .ckstop(f, "`x`, `y` and `z` must be numbers of the same length")
  }
  if (length(x) < 3L || !all(is.finite(c(x, y, z)))) {
    .ckstop(f, "`x`, `y` and `z` must be finite, at least three points")
  }
}
