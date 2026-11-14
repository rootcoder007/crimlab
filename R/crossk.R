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
#'   the Royal Statistical Society, Series B*, 39(2), 172-192.
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

#' Quadrat count test
#'
#' Counts the points in an `nx` by `ny` grid of equal tiles of the
#' rectangle `win` and tests complete spatial randomness by Pearson's
#' chi-square, `X2 = sum (n_k - e)^2 / e` with `e` the mean count, on `nx ny
#' - 1` degrees of freedom: the upper tail against clustering and the
#' two-sided p-value (twice the smaller tail, `quadrat.test()`'s default),
#' which also flags regular patterns. A point on an inner tile boundary
#' goes to the tile above or to its right, one on the outer edge to the
#' tile inside. These are the values of `spatstat.explore::quadrat.test()`.
#' The counts feed [morsta()].
#'
#' @inheritParams crossk
#' @param nx,ny The tiles across and up: whole numbers from 1, at least 2
#'   tiles.
#' @return A list of class `qdrtst`: `counts`, a matrix with one row per row
#'   of tiles from the top and one column per column from the left; `stat`,
#'   a named vector: `x2`, `df`, `p` (upper tail), `ptwo` and `lambda` (the
#'   points per unit area).
#' @references Greig-Smith, P. (1952). The use of random and contiguous
#'   quadrats in the study of the structure of plant communities. *Annals
#'   of Botany*, 16(2), 293-316.
#'   \doi{10.1093/oxfordjournals.aob.a083317}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(100, 0, 2), y = runif(100))
#' q <- qdrtst(p, win = c(0, 2, 0, 1), nx = 4, ny = 2)
#' q
#' morsta(q$counts)
#' @export
qdrtst <- function(x, y = NULL, win, nx, ny = nx) {
  .Call(C_clab_r_qdrtst, x, y, win, nx, ny)
}

#' @export
print.qdrtst <- function(x, ...) {
  s <- x[["stat"]]
  cat("qdrtst: X2 ", format(signif(s[["x2"]], 5)), " on ", s[["df"]],
      " df, p ", format(signif(s[["p"]], 3)), " (two-sided ",
      format(signif(s[["ptwo"]], 3)), ")\n", sep = "")
  print(x[["counts"]])
  invisible(x)
}

#' Morisita's index of dispersion
#'
#' Morisita's (1959) index of the counts in `n` equal quadrats, `I = n
#' (sum x^2 - S) / (S^2 - S)` with `S` the total: 1 for a random pattern,
#' above for clumping (`n` when all fall in one quadrat), below for an even
#' spread. The uniform and clumped indices are the values of `I` at the
#' `crit / 2` tails of chi-square on `n - 1` degrees of freedom, and the
#' standardised index of Smith-Gill (1975) maps `I` onto `[-1, 1]` with
#' them at `-0.5` and `0.5`; `chi = I (S - 1) + n - S` has the chi-square
#' distribution under randomness. These are the values of
#' `vegan::dispindmorisita()`.
#'
#' @param counts The counts: numbers from 0, at least two, adding up to
#'   more than 1 (a [qdrtst()] count matrix will do).
#' @param crit The level of the critical indices, above 0 and below 1.
#' @return A named vector: `imor`, `mclu`, `muni`, `imst`, `chi` and `p`
#'   (its upper chi-square tail).
#' @references Morisita, M. (1959). Measuring of the dispersion of
#'   individuals and analysis of the distributional patterns. *Memoirs of
#'   the Faculty of Science, Kyushu University, Series E (Biology)*, 2,
#'   215-235.
#'
#'   Smith-Gill, S. J. (1975). Cytophysiological basis of disruptive
#'   pigmentary patterns in the leopard frog *Rana pipiens*. II. Wild type
#'   and mutant cell-specific patterns. *Journal of Morphology*, 146(1),
#'   35-54. \doi{10.1002/jmor.1051460103}
#' @examples
#' morsta(c(3, 0, 7, 12, 1, 4, 0, 9))
#' @export
morsta <- function(counts, crit = 0.05) {
  .Call(C_clab_r_morsta, counts, crit)
}

#' Indices of residential segregation
#'
#' The evenness of groups across areal units (Massey and Denton 1988),
#' from a table of counts with one row per unit and one column per group:
#' Duncan's dissimilarity `D = sum |x_ik / X_k - x_il / X_l| / 2` between
#' each pair of groups (the share of either group that would have to move
#' for the two to be spread alike), and for each group against the rest
#' the segregation index `IS`, Theil's entropy index `H`, Atkinson's index
#' with shape `delta` and the Gini index. With weights `w` between the
#' units, Morrill's (1991) adjusted `D` subtracts the weighted mean
#' absolute difference of the pair's shares between neighbouring units, so
#' that segregation across a boundary counts less than segregation within.
#' Units with no one in them are left out. That `D` lies between 0 and 1 is
#' proved in Lean: see `proof("dissix")`.
#'
#' @param x The counts: a matrix or data frame of numbers from 0, units by
#'   groups (named columns), every group with a positive total.
#' @param delta Atkinson's shape, above 0 and below 1.
#' @param w `NULL`, or weights between the units (as [gmoran()] takes
#'   them), with a positive total.
#' @return A list of class `dissix`: `d`, the matrix of `D`; `grp`, a data
#'   frame of `group`, `is`, `h`, `atk` and `gini`; `dm`, `NULL` or the
#'   matrix of Morrill's `D`; and `delta`.
#' @references Duncan, O. D. and Duncan, B. (1955). A methodological
#'   analysis of segregation indexes. *American Sociological Review*,
#'   20(2), 210-217. \doi{10.2307/2088328}
#'
#'   Massey, D. S. and Denton, N. A. (1988). The dimensions of residential
#'   segregation. *Social Forces*, 67(2), 281-315. \doi{10.2307/2579183}
#'
#'   Morrill, R. L. (1991). On the measure of geographic segregation.
#'   *Geography Research Forum*, 11, 25-36.
#'   \doi{10.64166/p1b4ca43}
#' @examples
#' x <- cbind(a = c(90, 80, 10, 5, 50), b = c(10, 20, 85, 95, 50))
#' dissix(x)
#' dissix(x, w = spwt(nb = list(2, c(1, 3), c(2, 4), c(3, 5), 4),
#'                    type = "contig"))$dm
#' @export
dissix <- function(x, delta = 0.5, w = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_dissix, x, delta, w)
}

#' @export
print.dissix <- function(x, ...) {
  cat("dissix: dissimilarity between groups\n")
  print(x[["d"]])
  print(x[["grp"]], row.names = FALSE)
  if (!is.null(x[["dm"]])) {
    cat("Morrill's adjusted D\n")
    print(x[["dm"]])
  }
  invisible(x)
}
