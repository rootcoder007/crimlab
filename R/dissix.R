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
