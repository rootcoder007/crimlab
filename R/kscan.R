#' Space-time scan statistic of incidents
#'
#' Kulldorff's space-time permutation scan (Kulldorff, Heffernan, Hartman,
#' Assuncao and Mostashari 2005), which needs the incidents alone, no
#' population: cylinders are circles of radius `r` about each distinct
#' incident location times windows `[s, s + tw)` that start at each
#' distinct incident time. A cylinder holding `c` incidents, of `n_C` in its
#' circle (at any time) and `n_W` in its window (anywhere), expects `e =
#' n_C n_W / n` when place and time are independent, and scores the
#' Poisson log-likelihood ratio `c log(c / e) + (n - c) log((n - c) / (n -
#' e))` when `c > e`, 0 otherwise (Kulldorff 1997). The most likely
#' cluster has the largest ratio; secondary clusters follow in order of
#' their ratio, each the best of its circle among those whose circles
#' overlap none chosen before. Each p-value compares the cluster's ratio
#' with the largest ratio of the whole scan with the times shuffled among
#' the locations `nsim` times from `seed`, so a secondary cluster's p-value
#' is conservative. That the ratio is never negative is proved in Lean:
#' see `proof("kscan")`. The work grows as the number of locations times
#' `length(r)` times the distinct times times `length(tw)`, for the scan
#' and for each permutation.
#'
#' @inheritParams knxnr
#' @param r The radii of the circles: positive and increasing, in the unit
#'   of the coordinates (metres when `lonlat = TRUE`).
#' @param tw The lengths of the windows: positive, in the unit of `t`.
#' @param minc The fewest incidents a cluster holds: a whole number from 1.
#' @param top The most clusters reported: a whole number from 1.
#' @param nsim The number of permutations: a whole number from 0.
#' @return A list of class `kscan`: `cl`, a data frame with one row per
#'   cluster, most likely first: its centre `x`, `y` (in degrees when
#'   `lonlat = TRUE`), `r`, `start` and `end` (the window `[start, end)`),
#'   `n` (incidents in the cylinder), `exp`, `ncirc` and `nwin` (incidents
#'   in its circle and in its window), `rr` (the relative risk inside
#'   against outside), `llr` and `p` (`(1 + #{max* >= llr}) / (nsim + 1)`,
#'   `NaN` without permutations); `memb`, each incident's cluster or `NA`;
#'   `n` and `nsim`.
#' @references Kulldorff, M. (1997). A spatial scan statistic.
#'   *Communications in Statistics - Theory and Methods*, 26(6), 1481-1496.
#'   \doi{10.1080/03610929708831995}
#'
#'   Kulldorff, M., Heffernan, R., Hartman, J., Assuncao, R. and
#'   Mostashari, F. (2005). A space-time permutation scan statistic for
#'   disease outbreak detection. *PLoS Medicine*, 2(3), e59.
#'   \doi{10.1371/journal.pmed.0020059}
#' @examples
#' set.seed(1)
#' p <- cbind(x = c(runif(150, 0, 20), rnorm(20, 5, 0.3)),
#'            y = c(runif(150, 0, 20), rnorm(20, 15, 0.3)))
#' t <- c(sample(0:99, 150, TRUE), sample(40:44, 20, TRUE))
#' s <- kscan(p, t = t, r = c(1, 3), tw = c(5, 20), top = 2, nsim = 99,
#'            seed = 1)
#' s
#' @export
kscan <- function(x, y = NULL, t, r, tw, minc = 2, top = 1, nsim = 999,
                  seed = NULL, lonlat = FALSE) {
  if (inherits(t, "POSIXt")) {
    t <- as.numeric(as.POSIXct(t)) / 86400
  } else if (inherits(t, "Date")) {
    t <- as.numeric(t)
  }
  .Call(C_clab_r_kscan, x, y, t, r, tw, minc, top, nsim, seed, lonlat)
}

#' @export
print.kscan <- function(x, ...) {
  cat("kscan: ", nrow(x[["cl"]]), " clusters of ", x[["n"]],
      " incidents, ", x[["nsim"]], " permutations\n", sep = "")
  print(x[["cl"]], row.names = FALSE)
  invisible(x)
}
