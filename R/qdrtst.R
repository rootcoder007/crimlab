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
