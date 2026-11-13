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
