#' Cross-type K function and the segregation test
#'
#' Whether incidents of one type (`i`) have incidents of another (`j`)
#' near them more than independence gives: the cross-type K function
#' `K_ij(r) = A / (n_i n_j) sum e_ab 1(d_ab <= r)` over pairs of a point
#' `a` of type `i` and a point `b` of type `j` in the rectangle `win` of
#' area `A`, with Ripley's (1977) isotropic correction `e_ab` (the length of
#' the circle about `a` through `b` over its length inside `win`), and `L
#' = sqrt(K / pi)` (Lotwick and Silverman 1982). Under independence `K =
#' pi r^2`. With `i` equal to `j` it is Ripley's K of that type. These are
#' the values of `spatstat.explore::Kcross()` with `correction =
#' "isotropic"`.
#'
#' With `sigma`, the segregation test of Diggle, Zheng and Durr (2005)
#' across all the types: at each point the share of each type among the
#' other points, weighted by a normal kernel of standard deviation
#' `sigma`, against the type's overall share, `T = sum_points sum_types
#' (p_m - pbar_m)^2`, with a p-value from `nsim` random relabellings of
#' the types (`(1 + #{T* >= T}) / (nsim + 1)`), the statistic of
#' `spatstat.explore::segregation.test()`.
#'
#' @inheritParams kdsurf
#' @param marks The type of each point: labels with no missing values.
#' @param i,j The two types, by label or number (of the sorted labels, or
#'   a factor's levels).
#' @param win The study rectangle `c(xmin, xmax, ymin, ymax)`, holding
#'   every point.
#' @param r The distances: increasing, from 0.
#' @param sigma `NULL` for no segregation test, or the kernel's standard
#'   deviation, positive.
#' @param nsim The number of relabellings: a whole number from 0.
#' @inheritParams nnidx
#' @return A list of class `crossk`: `k`, a data frame of `r`, `k`, `l` and
#'   `theo` (`pi r^2`); `seg`, `NULL` or a named vector of `t`, `p` (`NaN`
#'   without relabellings) and `nsim`; `i` and `j`, the types' labels.
#' @references Ripley, B. D. (1977). Modelling spatial patterns. *Journal of
#'   the Royal Statistical Society, Series B*, 39(2), 172-212.
#'   \doi{10.1111/j.2517-6161.1977.tb01615.x}
#'
#'   Lotwick, H. W. and Silverman, B. W. (1982). Methods for analysing
#'   spatial processes of several types of points. *Journal of the Royal
#'   Statistical Society, Series B*, 44(3), 406-413.
#'   \doi{10.1111/j.2517-6161.1982.tb01221.x}
#'
#'   Diggle, P. J., Zheng, P. and Durr, P. (2005). Nonparametric estimation
#'   of spatial segregation in a multivariate point process: bovine
#'   tuberculosis in Cornwall, UK. *Journal of the Royal Statistical
#'   Society, Series C*, 54(3), 645-658.
#'   \doi{10.1111/j.1467-9876.2005.05373.x}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(120), y = runif(120))
#' m <- ifelse(p[, 1] < 0.5, "burglary", "robbery")
#' crossk(p, marks = m, i = "burglary", j = "robbery",
#'        win = c(0, 1, 0, 1), r = c(0.05, 0.1, 0.2), sigma = 0.1,
#'        nsim = 19, seed = 1)
#' @export
crossk <- function(x, y = NULL, marks, i, j, win, r, sigma = NULL, nsim = 99,
                   seed = NULL) {
  .Call(C_clab_r_crossk, x, y, marks, i, j, win, r, sigma, nsim, seed)
}

#' @export
print.crossk <- function(x, ...) {
  cat("crossk: K of ", x[["j"]], " about ", x[["i"]], "\n", sep = "")
  print(x[["k"]], row.names = FALSE)
  s <- x[["seg"]]
  if (!is.null(s)) {
    cat("segregation T ", format(signif(s[["t"]], 5)), ", p ",
        format(signif(s[["p"]], 3)), " (", s[["nsim"]], " relabellings)\n",
        sep = "")
  }
  invisible(x)
}
