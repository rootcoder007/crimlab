#' Spatial weights
#'
#' Builds the neighbours and weights the autocorrelation statistics
#' ([gearc()], [gorg()], [lmorn()], [locgi()]) use: the `k` nearest
#' neighbours of each point (`type = "knn"`, ties broken by the lower row
#' number), the points within distance `d` (`"band"`, the distance itself
#' included), or a contiguity list `nb` (`"contig"`, one element per unit
#' giving its neighbours' numbers, as spdep's `poly2nb()` returns). Weights
#' are 1, or the inverse distance to the power `idw`; `style = "row"`
#' divides each unit's weights by their sum (a unit without neighbours keeps
#' none), `"raw"` leaves them.
#'
#' @inheritParams cgstat
#' @param type `"knn"`, `"band"` or `"contig"`.
#' @param k For `"knn"`: the number of neighbours, from 1 to the number of
#'   points less 1.
#' @param d For `"band"`: the distance, positive, in the unit of the
#'   coordinates (metres when `lonlat = TRUE`).
#' @param nb For `"contig"`: a list with one element per unit, each a vector
#'   of neighbour numbers (from 1; `0` alone for none).
#' @param idw The power of the inverse distance: 0 (binary weights) or more;
#'   above 0 it needs the points.
#' @param style `"row"` or `"raw"`.
#' @param metric `"euclid"` or `"manhattan"` (`|dx| + |dy|`).
#' @return A list: `n` (the number of units), and one element per
#'   neighbour pair, in order of unit then neighbour: `i`, `j` (integer unit
#'   numbers), `d` (their distance, `NA` for contiguity without points) and
#'   `w` (the weight). `x` may be `NULL` for contiguity with binary weights.
#' @references Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes:
#'   Models and Applications*. Pion, London.
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(20), y = runif(20))
#' w <- spwt(p, k = 3)
#' str(w)
#' spwt(nb = list(2, c(1, 3), 2), type = "contig")$w
#' @export
spwt <- function(x = NULL, y = NULL, type = "knn", k = 4, d = NULL,
                 nb = NULL, idw = 0, style = "row", metric = "euclid",
                 lonlat = FALSE) {
  .Call(C_clab_r_spwt, x, y, type, k, d, nb, idw, style, metric, lonlat)
}

#' Geary's C
#'
#' Geary's (1954) contiguity ratio of values on units with spatial weights
#' `w`: `C = (n - 1) sum w_ij (x_i - x_j)^2 / (2 S0 sum (x_i - xbar)^2)`, 1
#' under no autocorrelation, below 1 when neighbours are alike. Its variance
#' is given under normality and under randomisation (Cliff and Ord 1981),
#' with z statistics `(1 - C) / sd` and their upper normal tails, and, with
#' `nsim > 0`, a permutation p-value against positive autocorrelation.
#' That C has expectation 1 under randomisation is proved in Lean: see
#' `proof("gearc")`.
#'
#' @param x The values: finite numbers, one per unit.
#' @param w The weights: [spwt()]'s list, or an spdep `listw` object, for as
#'   many units as `x`.
#' @inheritParams gstr
#' @inheritParams nnidx
#' @return A named numeric vector: `c`, `vn` and `vr` (its variance under
#'   normality and randomisation), `zn`, `zr`, `pn`, `pr`, `pperm`
#'   (`(1 + #{C* <= C}) / (nsim + 1)`, `NaN` without permutations) and `s0`
#'   (the sum of the weights).
#' @references Geary, R. C. (1954). The contiguity ratio and statistical
#'   mapping. *The Incorporated Statistician*, 5(3), 115-145.
#'   \doi{10.2307/2986645}
#'
#'   Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
#'   Applications*. Pion, London.
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(50), y = runif(50))
#' gearc(p[, 1] + p[, 2] + rnorm(50, sd = 0.1), spwt(p, k = 4))
#' @export
gearc <- function(x, w, nsim = 0, seed = NULL) {
  .Call(C_clab_r_gearc, x, w, nsim, seed)
}

#' The general G of Getis and Ord
#'
#' The global concentration statistic of Getis and Ord (1992) for
#' non-negative values: `G = sum_ij w_ij x_i x_j / sum_{i != j} x_i x_j`,
#' its expectation `S0 / (n (n - 1))` and variance under randomisation, the
#' z statistic, its upper normal tail (high values near high values) and,
#' with `nsim > 0`, a permutation p-value. Unlike Moran's I and Geary's C it
#' tells high-value clusters from low-value ones.
#' Its expectation under randomisation is proved in Lean: see
#' `proof("gorg")`.
#'
#' @inheritParams gearc
#' @param x The values: finite numbers from 0, at least two positive.
#' @return A named numeric vector: `g`, `e`, `v`, `z`, `p` and `pperm`
#'   (`(1 + #{G* >= G}) / (nsim + 1)`).
#' @references Getis, A. and Ord, J. K. (1992). The analysis of spatial
#'   association by use of distance statistics. *Geographical Analysis*,
#'   24(3), 189-206. \doi{10.1111/j.1538-4632.1992.tb00261.x}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(50), y = runif(50))
#' gorg(exp(3 * p[, 1]), spwt(p, type = "band", d = 0.2, style = "raw"))
#' @export
gorg <- function(x, w, nsim = 0, seed = NULL) {
  .Call(C_clab_r_gorg, x, w, nsim, seed)
}

#' Local Moran (LISA)
#'
#' Anselin's (1995) local indicator of spatial association for each unit:
#' with `z = x - mean(x)` and `m2 = mean(z^2)`, `I_i = z_i / m2 sum_j w_ij
#' z_j`, its expectation and variance under total randomisation, the
#' two-sided normal p-value, and a conditional permutation p-value: unit `i`
#' keeps its value and its neighbours' values are drawn without replacement
#' from the others, `nsim` times from `seed`. The quadrant tells
#' high-high and low-low clusters from high-low and low-high outliers. The
#' p-values are per unit and not adjusted for the many tests made (use
#' [stats::p.adjust()]); with row-standardised weights the sum of the `I_i`
#' is `S0` times global Moran's I. That identity is proved in Lean: see
#' `proof("lmorn")`.
#'
#' @inheritParams gearc
#' @return A data frame with one row per unit: `ii` (I_i), `e`, `v`, `z`,
#'   `p`, `pperm` (`(1 + min(#{I* >= I}, #{I* <= I})) / (nsim + 1)`, `NaN`
#'   without permutations) and `quad` (1 high-high, 2 low-low, 3 low-high,
#'   4 high-low, 0 neither).
#' @references Anselin, L. (1995). Local indicators of spatial association:
#'   LISA. *Geographical Analysis*, 27(2), 93-115.
#'   \doi{10.1111/j.1538-4632.1995.tb00338.x}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(60), y = runif(60))
#' v <- ifelse(p[, 1] > 0.7 & p[, 2] > 0.7, 5, 1) + rnorm(60, sd = 0.3)
#' l <- lmorn(v, spwt(p, k = 5), nsim = 199, seed = 1)
#' head(l)
#' @export
lmorn <- function(x, w, nsim = 0, seed = NULL) {
  .Call(C_clab_r_lmorn, x, w, nsim, seed)
}

#' Local G and G* of Getis and Ord
#'
#' The local statistics of Getis and Ord (1992) as z values (Ord and Getis
#' 1995): `G_i` compares the weighted sum of the neighbours' values with its
#' expectation from the other `n - 1` values; `G*_i` (`star = TRUE`)
#' includes the unit itself with weight 1 and uses all `n`. Large positive
#' values mark concentrations of high values, large negative ones of low
#' values. With `star = TRUE` the weights are best binary (`style = "raw"`).
#'
#' @inheritParams gearc
#' @param star `TRUE` for G*, `FALSE` for G.
#' @return A numeric vector of z values, one per unit; `NaN` for a unit
#'   without neighbours.
#' @references Ord, J. K. and Getis, A. (1995). Local spatial
#'   autocorrelation statistics: distributional issues and an application.
#'   *Geographical Analysis*, 27(4), 286-306.
#'   \doi{10.1111/j.1538-4632.1995.tb00912.x}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(60), y = runif(60))
#' v <- rpois(60, ifelse(p[, 1] > 0.7, 6, 2))
#' locgi(v, spwt(p, type = "band", d = 0.2, style = "raw"))
#' @export
locgi <- function(x, w, star = TRUE) {
  .Call(C_clab_r_locgi, x, w, star)
}

#' Directional mean and circular variance
#'
#' The mean direction of a set of directions, the mean resultant length
#' `R` (from 0, directions spread evenly, to 1, all the same), the circular
#' variance `1 - R` and the circular standard deviation
#' `sqrt(-2 log R)` (Mardia and Jupp 2000). Directions are given as angles,
#' or as points seen from an origin (the directions of incidents from an
#' offender's base, say). The mean direction of points depends on the
#' origin chosen. That `R` lies in `[0, 1]` is proved in Lean: see
#' `proof("dirmn")`.
#'
#' @inheritParams cgstat
#' @param origin For points: the origin `c(x, y)` (in degrees when
#'   `lonlat = TRUE`); no point may lie on it.
#' @param w `NULL`, or one weight from 0 per direction, positive sum.
#' @param angle Instead of points: angles in radians, anticlockwise from
#'   east (the x axis).
#' @return A named numeric vector: `theta` (the mean direction in radians
#'   anticlockwise from east, in `[0, 2 pi)`; `NaN` when `R` is 0),
#'   `bearing` (the same as a compass bearing, degrees clockwise from
#'   north), `r`, `var`, `sd`, `w` (the total weight) and `n`.
#' @references Mardia, K. V. and Jupp, P. E. (2000). *Directional
#'   Statistics*. Wiley. \doi{10.1002/9780470316979}
#' @examples
#' dirmn(c(1, 2, 0), c(1, 1, 3), origin = c(0, 0))
#' dirmn(angle = c(0.1, 0.2, 6.2))
#' @export
dirmn <- function(x = NULL, y = NULL, origin = c(0, 0), w = NULL,
                  angle = NULL, lonlat = FALSE) {
  .Call(C_clab_r_dirmn, x, y, origin, w, angle, lonlat)
}
