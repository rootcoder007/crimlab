#' Difference in differences with two-way fixed effects
#'
#' The effect of a policy on a unit-by-period panel (towns by year, beats by
#' month, courts by quarter) from a regression of the outcome on the
#' treatment with unit and period fixed effects, with standard errors
#' clustered by unit or by any coarser grouping, one-way or two-way. Counts
#' can be fitted by Poisson or negative binomial regression with the same
#' effects. Options give the estimate leaving out each treated unit in
#' turn, its Goodman-Bacon decomposition, and the estimate after matching
#' treated units to comparison units on their outcomes before treatment.
#'
#' The model is `y_it = b trt_it + adj_it'g + a_i + m_t + e_it`, fitted by
#' least squares after the within transformation (alternating projections,
#' so an unbalanced panel is fine; weights `w` give weighted least squares
#' with weighted means). With `fe = "time"` the unit effects are left out.
#' The variance is the cluster-robust CR1 estimator, `G / (G - 1) (n - 1) /
#' (n - K)` times the sandwich, `K` counting the regressors, the constant,
#' and the levels less one of each fixed effect not nested in a cluster;
#' with two kinds of cluster it is `V1 + V2 - V12` (Cameron, Gelbach and
#' Miller 2011) with `G` the smaller number of clusters. Tests and intervals
#' use Student's t on `G - 1` degrees of freedom. These are the conventions
#' of 'fixest' and Stata's 'reghdfe'. Fewer than ten clusters give a
#' warning: the wild cluster bootstrap is the tool for so few.
#'
#' With `fam = "pois"` or `"nb"` the model is `E(y_it) = exp(b trt_it +
#' adj_it'g + a_i + m_t)` times the exposure `expo`, fitted by maximum
#' likelihood with dummies for the effects (the steps of [cntreg()]). Units
#' and periods whose counts are all zero carry nothing about `b` and are
#' dropped first. The Poisson estimate is consistent whatever the
#' dispersion (Wooldridge 1999); the sandwich and its factor are those of
#' the least-squares fit, and intervals are normal; `rr` is the rate ratio
#' `exp(b)`. With unit dummies the negative binomial is the unconditional
#' fixed-effects model, which behaves well when there are many periods per
#' unit (Allison and Waterman 2002).
#'
#' `releff` is the estimate over `ctlmean`, the mean outcome of the units
#' never treated before the first treatment (the untreated rows before it
#' when every unit is treated at some point): "a 12% fall". For counts it is
#' `rr - 1`.
#'
#' With `bacon = TRUE` (a balanced panel, treatment that stays on, no
#' adjustment or weights), the estimate is split into the two-group,
#' two-period comparisons it averages, with their weights (Goodman-Bacon
#' 2021, theorem 1): treated against never treated, earlier against later
#' treated before the later start, and later against earlier treated after
#' the earlier start (already-treated units as controls). The weights sum
#' to 1 and `sum(weight * est)` is the estimate exactly. Units treated from
#' the first period act as an already-treated group.
#'
#' With `mtch` a whole number above 0, each treated unit is matched, with
#' replacement, to its `mtch` nearest units never treated by the
#' Mahalanobis distance of their outcomes in every period before the first
#' adoption, with the pooled within-group covariance (as 'MatchIt' computes
#' it), within its `exact` group and, with `calip`, within `calip` standard
#' deviations of all the distances computed. A comparison unit's weight is
#' the sum, over the treated units it serves, of one over their number of
#' matches; matched treated units weigh 1 and the rest 0. The estimate is
#' then the weighted fit. Matching on the outcomes before treatment can
#' bring in regression to the mean when units were chosen for an unusual
#' pre-period (Daw and Hatfield 2018).
#'
#' @param y The outcome: finite numbers, one per row (counts from 0 for
#'   `fam = "pois"` or `"nb"`).
#' @param unit The unit of each row: labels.
#' @param time The period of each row: numbers (a year, a month index);
#'   each unit at most once per period.
#' @param trt The treatment of each row: 0 or 1 (or `FALSE` or `TRUE`).
#' @param adj `NULL` or time-varying controls: a data frame or list of
#'   columns, a numeric matrix, or one vector, one value per row. Factor,
#'   text and logical columns enter as indicators of each level but the
#'   first.
#' @param clus `NULL` (cluster by unit), the cluster of each row (labels),
#'   or a list or data frame of two such columns for two-way clustering.
#' @param w `NULL` or weights from 0, one per row (population, matching
#'   weights); rows of weight 0 take no part.
#' @param fe `"both"` (unit and period effects) or `"time"` (period effects
#'   only).
#' @param fam `"ols"` (least squares), `"pois"` (Poisson) or `"nb"`
#'   (negative binomial) regression.
#' @param expo `NULL` or the exposure of each row for counts (population,
#'   person-days): positive numbers; its logarithm is the offset.
#' @param loo `TRUE` for the estimate leaving out each ever-treated unit in
#'   turn.
#' @param bacon `TRUE` for the Goodman-Bacon decomposition.
#' @param mtch `0`, or the number of comparison units matched to each
#'   treated unit before the fit.
#' @param exact `NULL` or labels, constant within units: match only within
#'   them (a region, a district).
#' @param calip `NULL` or the largest distance of a match, in standard
#'   deviations of the distances.
#' @param level The confidence level of the intervals: one number strictly
#'   between 0 and 1.
#' @return An object of class `didfe`, a list: `coef`, one row per term
#'   (the treatment `"trt"` first, then the controls): `term`, `est`, `se`,
#'   `lo`, `hi`, `stat` (t or z), `p`, and for counts `rr`, `rrlo`, `rrhi`;
#'   `fit`, one row: `n` (rows used), `units`, `periods`, `ncl` (`G`), `k`
#'   (`K`), `ctlmean`, `releff`, `theta` (negative binomial), `dropped`
#'   (units with no counts); `vcov`, the variance of the terms; `loo`
#'   (`unit`, `est`, `se`), `bacon` (`type`, `treated`, `control`, `est`,
#'   `weight`; `control` is `Inf` for the never treated) and `match`
#'   (`treated`, `control`, `dist`), each `NULL` unless asked for; `fam`;
#'   and `level`.
#' @references
#'   Bertrand, M., Duflo, E. and Mullainathan, S. (2004). How much should
#'   we trust differences-in-differences estimates? *Quarterly Journal of
#'   Economics*, 119(1), 249-275. \doi{10.1162/003355304772839588}
#'
#'   Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2011). Robust
#'   inference with multiway clustering. *Journal of Business and Economic
#'   Statistics*, 29(2), 238-249. \doi{10.1198/jbes.2010.07136}
#'
#'   Goodman-Bacon, A. (2021). Difference-in-differences with variation in
#'   treatment timing. *Journal of Econometrics*, 225(2), 254-277.
#'   \doi{10.1016/j.jeconom.2021.03.014}
#'
#'   Wooldridge, J. M. (1999). Distribution-free estimation of some
#'   nonlinear panel data models. *Journal of Econometrics*, 90(1), 77-97.
#'   \doi{10.1016/S0304-4076(98)00033-5}
#'
#'   Allison, P. D. and Waterman, R. P. (2002). Fixed-effects negative
#'   binomial regression models. *Sociological Methodology*, 32, 247-265.
#'   \doi{10.1111/1467-9531.00117}
#'
#'   Daw, J. R. and Hatfield, L. A. (2018). Matching and regression to the
#'   mean in difference-in-differences analysis. *Health Services
#'   Research*, 53(6), 4138-4156. \doi{10.1111/1475-6773.12993}
#' @examples
#' # twelve towns over ten years; four adopt a by-law in 2019, four in 2021
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
#' i <- match(d$town, LETTERS)
#' first <- c(2019, 2019, 2019, 2019, 2021, 2021, 2021, 2021, NA, NA, NA, NA)
#' d$bylaw <- as.integer(!is.na(first[i]) & d$year >= first[i])
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * d$bylaw + sin(seq_len(120))
#' r <- didfe(d$rate, d$town, d$year, d$bylaw)
#' r
#' # its comparisons, and each adopting town left out in turn
#' didfe(d$rate, d$town, d$year, d$bylaw, bacon = TRUE, loo = TRUE)$bacon
#' @export
didfe <- function(y, unit, time, trt, adj = NULL, clus = NULL, w = NULL,
                  fe = "both", fam = "ols", expo = NULL, loo = FALSE,
                  bacon = FALSE, mtch = 0, exact = NULL, calip = NULL,
                  level = 0.95) {
  .Call(C_clab_r_didfe, y, unit, time, trt, adj, clus, w, fe, fam, expo, loo,
        bacon, mtch, exact, calip, level)
}

#' @export
print.didfe <- function(x, ...) {
  f <- x[["fit"]]
  cat("didfe: ", x[["fam"]], ", ", f[["n"]], " rows, ", f[["units"]],
      " units, ", f[["periods"]], " periods, ", f[["ncl"]], " clusters\n",
      sep = "")
  print(x[["coef"]], row.names = FALSE)
  cat("comparison mean ", format(f[["ctlmean"]]), ", relative effect ",
      format(f[["releff"]]), "\n", sep = "")
  if (!is.null(x[["bacon"]])) {
    b <- x[["bacon"]]
    cat("Goodman-Bacon weights by type:\n")
    print(tapply(b[["weight"]], b[["type"]], sum))
  }
  invisible(x)
}
