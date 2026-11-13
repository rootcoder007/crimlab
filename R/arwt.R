#' Crosswalk of counts between two sets of areal units
#'
#' Moves values (counts, populations) kept on one set of units (an old
#' boundary version, say) onto another, through the pairs of a source and
#' a target unit that overlap: each source's value is split among its
#' targets in proportion to `area`, the size of the overlap (its area, or
#' its population for a population-weighted crosswalk), which assumes the
#' value spread uniformly over the source. Where each source lies in one
#' target (a boundary version that only merges units) the result is exact
#' aggregation and `exact` is `TRUE`; splitting a source is an
#' approximation. A source with no pair is lost (`lost`), and each source's
#' main target (its largest share) maps an old identifier onto the new
#' scheme. The concentration of the result across the targets is reported
#' as the Gini coefficient and the share of the largest `top` targets.
#'
#' @param x The values of the source units: finite numbers.
#' @param from,to The overlapping pairs: the source unit (a label in `id`)
#'   and the target unit of each, as many as there are pairs.
#' @param area `NULL` to split each source evenly among its targets, or the
#'   size of each overlap, from 0.
#' @param id The labels of the source units, one per value of `x`, each
#'   once; by default `names(x)`.
#' @param top The numbers of largest targets whose share is reported:
#'   whole numbers from 1.
#' @return A list of class `arwt`: `out`, a data frame of `to` and `value`
#'   (every target in a pair, in sorted order); `map`, a data frame of
#'   `from`, `to` (its main target) and `share`; `exact`; `lost`; and
#'   `conc`, a named vector: `gini` and `top<k>`.
#' @references Goodchild, M. F. and Lam, N. S.-N. (1980). Areal
#'   interpolation: a variant of the traditional spatial problem.
#'   *Geo-Processing*, 1, 297-312.
#' @examples
#' # four old units onto three new ones; unit d is split 70/30 by area
#' arwt(c(a = 120, b = 40, c = 60, d = 100),
#'      from = c("a", "b", "c", "d", "d"), to = c("N1", "N1", "N2", "N2", "N3"),
#'      area = c(1, 1, 1, 7, 3), top = 1)
#' @export
arwt <- function(x, from, to, area = NULL, id = names(x), top = c(10, 20)) {
  if (is.null(id)) {
    stop("arwt: `id` must name each source unit (or `x` have names).",
         call. = FALSE)
  }
  .Call(C_clab_r_arwt, x, id, from, to, area, top)
}

#' @export
print.arwt <- function(x, ...) {
  cat("arwt: ", nrow(x[["map"]]), " sources onto ", nrow(x[["out"]]),
      " targets, ", if (x[["exact"]]) "exact (no source split)" else
        "sources split by overlap", if (x[["lost"]] != 0)
          paste0(", ", format(x[["lost"]]), " lost"), "\n", sep = "")
  print(x[["out"]], row.names = FALSE)
  print(signif(x[["conc"]], 4))
  invisible(x)
}
