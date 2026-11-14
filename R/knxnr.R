#' Knox test of near repeats
#'
#' Asks whether incidents close together in space also happen close
#' together in time more often than chance gives: the near-repeat pattern
#' of burglary and shootings (Knox 1964; Townsley, Homel and Chaseling
#' 2003). Every pair of incidents within the largest distance is put in a
#' band of distance, `(ds[k - 1], ds[k]]`, and a band of time, `(dt[l - 1],
#' dt[l]]` (the first bands start at 0 and include it), and each cell's
#' count is compared with its expectation when the times are exchangeable
#' among the locations, which is exact: `E = S_k T_l / P`, with `S_k` the
#' pairs in distance band `k`, `T_l` the pairs in time band `l` and `P =
#' n (n - 1) / 2`. Each cell has a one-sided permutation p-value, the times
#' shuffled among the locations `nsim` times from `seed`, and the Poisson
#' tail `P(X >= obs)` with the expectation as its mean, the classical Knox
#' approximation for a rare cell. Cells are tested separately: a table of
#' many cells shows some small p-values by chance alone. A near-repeat
#' pattern also arises when one offender commits a series; the test
#' speaks to timing, not to its cause.
#'
#' @inheritParams cgstat
#' @param t The times: numbers (days, say), or `Date` or `POSIXct` (taken
#'   in days), one per incident.
#' @param ds The upper edges of the distance bands: positive and
#'   increasing, in the unit of the coordinates (metres when `lonlat =
#'   TRUE`).
#' @param dt The upper edges of the time bands: positive and increasing, in
#'   the unit of `t`.
#' @param nsim The number of permutations: a whole number from 0 (none).
#' @inheritParams nnidx
#' @return A list of class `knxnr`: `tab`, a data frame with one row per
#'   cell, distance bands varying fastest: `dlo`, `dhi`, `tlo`, `thi`,
#'   `obs` (pairs), `exp`, `ratio` (`obs / exp`, `NaN` when nothing is
#'   expected), `p` (`(1 + #{obs* >= obs}) / (nsim + 1)`, `NaN` without
#'   permutations) and `ppois`; `n` (incidents), `pairs` (pairs within the
#'   largest distance) and `nsim`.
#' @references Knox, E. G. (1964). The detection of space-time
#'   interactions. *Journal of the Royal Statistical Society, Series C*,
#'   13(1), 25-30. \doi{10.2307/2985220}
#'
#'   Townsley, M., Homel, R. and Chaseling, J. (2003). Infectious
#'   burglaries: a test of the near repeat hypothesis. *British Journal of
#'   Criminology*, 43(3), 615-633.
#'   \doi{10.1093/bjc/43.3.615}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(150, 0, 2000), y = runif(150, 0, 2000))
#' k <- knxnr(p, t = runif(150, 0, 365), ds = c(200, 400), dt = c(14, 28),
#'            nsim = 99, seed = 2)
#' k
#' @export
knxnr <- function(x, y = NULL, t, ds, dt, nsim = 999, seed = NULL,
                  lonlat = FALSE) {
  if (inherits(t, "POSIXt")) {
    t <- as.numeric(as.POSIXct(t)) / 86400
  } else if (inherits(t, "Date")) {
    t <- as.numeric(t)
  }
  .Call(C_clab_r_knxnr, x, y, t, ds, dt, nsim, seed, lonlat)
}

#' @export
print.knxnr <- function(x, ...) {
  cat("knxnr: ", x[["n"]], " incidents, ", format(x[["pairs"]]),
      " close pairs, ", x[["nsim"]], " permutations\n", sep = "")
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}
