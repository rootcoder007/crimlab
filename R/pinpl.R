#' Points in polygons, and the population the units cover
#'
#' Finds the region holding each point (a facility, a station, an
#' incident) by the even-odd rule, a point on a boundary counting as
#' inside; a region given as several rings (parts and holes) is one
#' region, and a hole is any ring inside another of the same region. A
#' point outside every region gets `NA`; one inside overlapping regions
#' takes the first by label and `nreg` says how many hold it. Coordinates
#' are taken as planar (longitude and latitude work for regions small
#' enough that their edges are straight in degrees).
#'
#' With the regions' populations, the coverage: the population of the
#' regions holding a unit, and its share of the total. That share is not a
#' rate denominator. A region without a unit is not unserved: units serve
#' catchments, which a location does not state, and a rate built on the
#' regions that hold units pairs a partial denominator with a whole
#' numerator, inflating it unevenly.
#'
#' Polygons from a spatial file convert with, for instance,
#' `sf::st_coordinates()`, whose columns `X`, `Y`, `L1` (ring) and `L2` or
#' `L3` (feature) become `x`, `y`, `ring` and `id`.
#'
#' @inheritParams kdsurf
#' @param poly The regions: a data frame with columns `x` and `y` (or `lon`
#'   and `lat`), the vertices in order around each ring, `id` (the region's
#'   label) and, for regions of several rings, `ring`.
#' @param pop `NULL`, or the population of each region (from 0), in the
#'   sorted order of the labels.
#' @return A list of class `pinpl`: `reg` (each point's region label, `NA`
#'   outside), `nreg` (the regions holding it), `regs` (a data frame of
#'   `region`, `units`, `pop`, `hasunit` and `popshr`), and `cov`, a named
#'   vector: `regions`, `withunit`, `units`, `totpop`, `covpop`, `uncpop` and
#'   `covshr` (`NA` without populations).
#' @references Shimrat, M. (1962). Algorithm 112: position of point
#'   relative to polygon. *Communications of the ACM*, 5(8), 434.
#'   \doi{10.1145/368637.368653}
#' @examples
#' sq <- function(x0, id) {
#'   data.frame(x = x0 + c(0, 1, 1, 0), y = c(0, 0, 1, 1), id = id)
#' }
#' poly <- rbind(sq(0, "A"), sq(1, "B"), sq(2, "C"))
#' u <- pinpl(cbind(x = c(0.5, 0.2, 2.5, 9), y = c(0.5, 0.7, 0.5, 9)),
#'            poly = poly, pop = c(1200, 800, 450))
#' u
#' u$reg
#' @export
pinpl <- function(x, y = NULL, poly, pop = NULL) {
  .Call(C_clab_r_pinpl, x, y, poly, pop)
}

#' @export
print.pinpl <- function(x, ...) {
  v <- x[["cov"]]
  cat("pinpl: ", v[["units"]], " of ", length(x[["reg"]]),
      " points in ", v[["withunit"]], " of ", v[["regions"]], " regions",
      sep = "")
  if (!is.na(v[["covshr"]])) {
    cat("; they cover ", format(signif(100 * v[["covshr"]], 3)),
        "% of the population, which is not a rate denominator", sep = "")
  }
  cat("\n")
  print(x[["regs"]], row.names = FALSE)
  invisible(x)
}

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
#' @references Lam, N. S.-N. (1983). Spatial interpolation methods: a
#'   review. *The American Cartographer*, 10(2), 129-150.
#'   \doi{10.1559/152304083783914958}
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
