#' Treated, comparison and excluded places around sites
#'
#' The buffer design of place-based evaluations: a place (a block, an
#' address) within distance `r` of a site (a new facility, a treated
#' corner) is treated; one within `r` of a comparison point and farther
#' than `buf` from every site is a comparison; the rest are excluded, so
#' that spillover from the sites (displacement or a diffusion of benefits,
#' Bowers and Johnson 2003) does not reach the comparison. `buf` should be
#' several times `r`.
#'
#' @inheritParams cgstat
#' @param site,comp The sites and the comparison points, as coordinates in
#'   the form of `x` (longitude and latitude when `lonlat = TRUE`).
#' @param r The treatment radius: positive, in the unit of the coordinates
#'   (metres when `lonlat = TRUE`).
#' @param buf The distance from every site a comparison place keeps, from
#'   0.
#' @return A list of class `bufas`: `pts`, a data frame with one row per
#'   place: `grp` (a factor: treated, comparison, excluded), `dsite` and
#'   `dcomp` (the distances to the nearest site and comparison point); `n`,
#'   the counts of the groups; `par`, `r` and `buf`.
#' @references Bowers, K. J. and Johnson, S. D. (2003). Measuring the
#'   geographical displacement and diffusion of benefit effects of crime
#'   prevention activity. *Journal of Quantitative Criminology*, 19(3),
#'   275-301. \doi{10.1023/A:1024909009240}
#' @examples
#' set.seed(4)
#' p <- cbind(x = runif(200, 0, 5000), y = runif(200, 0, 5000))
#' b <- bufas(p, site = cbind(c(1000, 3500), c(1200, 3000)),
#'            comp = cbind(c(2500, 4500), c(500, 4000)), r = 330)
#' b
#' @export
bufas <- function(x, y = NULL, site, comp, r, buf = 3 * r, lonlat = FALSE) {
  .Call(C_clab_r_bufas, x, y, site, comp, r, buf, lonlat)
}

#' @export
print.bufas <- function(x, ...) {
  cat("bufas: radius ", format(x[["par"]][["r"]]), ", buffer ",
      format(x[["par"]][["buf"]]), "\n", sep = "")
  print(x[["n"]])
  invisible(x)
}
