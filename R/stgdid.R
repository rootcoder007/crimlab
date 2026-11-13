#' Effects of a policy adopted at different times
#'
#' Difference-in-differences estimators for staggered adoption that never
#' use an already-treated unit as a control: the group-time effects of
#' Callaway and Sant'Anna (2021), optionally with covariates; the
#' regression on cohort-period cells of Wooldridge (2021); and, for a
#' treatment that switches on and off, the estimator of de Chaisemartin and
#' D'Haultfoeuille (2020) with its placebo. Standard errors come from the
#' cluster bootstrap.
#'
#' Let cohort `g` be the units first treated in period `g`, and the
#' comparison units those never treated (`ctrl = "never"`) or not yet
#' treated by period `max(t, b)` (`"notyet"`, outside cohort `g`). For
#' `meth = "csa"`, `ATT(g, t)` is the cohort's mean change in the outcome
#' from the base period `b` less the comparison units', with `b = g - 1`
#' from adoption on and `b = t - 1` before (the `"varying"` base of the
#' 'did' package), so an effect before adoption is a placebo over one
#' period. With `adj`, parallel trends need hold only given the covariates,
#' taken at the base period, and each cell is the estimator of Sant'Anna
#' and Zhao (2020) on the cohort and its comparison units: `"dr"`, the
#' doubly robust `mean_D(dy - m(x)) - sum w (dy - m(x)) / sum w` with `m`
#' the least-squares fit of the change among the comparison units and `w =
#' p / (1 - p)` from the logistic propensity `p` of belonging to the
#' cohort; `"ipw"` without `m`; `"reg"` without the weights. Comparison
#' units with `p` of 0.995 or more are left out, as 'DRDID' does. Use
#' covariates fixed before treatment.
#'
#' For `meth = "wool"`, the regression of the outcome on unit and period
#' effects and one indicator for every treated cohort-period cell gives
#' each `ATT(g, t)` after adoption as a coefficient (without covariates
#' this is the imputation estimator); without units never treated, the
#' periods from the last adoption on are left out.
#'
#' The overall effect averages the `ATT(g, t)` after adoption weighted by
#' cohort size (the `"simple"` aggregation of 'did'); the effect at event
#' time `e` averages `ATT(g, g + e)` over the cohorts observed at `e`,
#' weighted by cohort size. Event time counts observed periods. Units
#' treated from the first period are left out.
#'
#' For `meth = "dcdh"`, in each period the units that switch into
#' treatment are compared with the units untreated in both periods, and
#' the units that switch out with the units treated in both: `DID_M =
#' sum_t (N_join,t DID_join,t + N_leave,t DID_leave,t) / N_S` over the
#' `N_S` switchers that have a comparison. The placebo runs the same
#' comparison one period earlier, for units whose treatment did not change
#' over the two periods before; it should be near zero.
#'
#' Standard errors are the standard deviations over `reps` bootstrap
#' samples of the clusters (units by default), drawn with replacement, with
#' everything recomputed; intervals are normal.
#'
#' @inheritParams didfe
#' @param first The period each row's unit is first treated (constant
#'   within units; `NA`, 0 or `Inf` for never), or `NULL` with `trt`.
#' @param trt The treatment of each row, 0 or 1, or `NULL` with `first`.
#'   For `meth = "csa"` and `"wool"` it must stay on once on.
#' @param meth `"csa"`, `"wool"` or `"dcdh"`.
#' @param ctrl `"notyet"` or `"never"`: the comparison units of `"csa"`.
#' @param adj `NULL` or covariates for `"csa"`: a data frame or list of
#'   columns, a numeric matrix or one vector, one value per row.
#' @param est With `adj`: `"dr"`, `"ipw"` or `"reg"`.
#' @param clus `NULL` (resample units) or the cluster of each row, constant
#'   within units.
#' @param reps The bootstrap samples: a whole number from 0 (no standard
#'   errors).
#' @param seed The seed of the bootstrap: `NULL` to draw it from R's
#'   generator (so [set.seed()] fixes it), or a whole number from 0.
#' @return An object of class `stgdid`, a list: `att`, the group-time
#'   effects (`cohort`, `period`, `evtm`, `est`, `se`, `lo`, `hi`, `ncoh`
#'   the cohort's size; `NULL` for `"dcdh"`); `overall` (`term`, `est`,
#'   `se`, `lo`, `hi`, and for `"dcdh"` the rows `"effect"` and `"placebo"`
#'   with `n`, the switchers); `dynamic`, by event time (`evtm`, `est`,
#'   `se`, `lo`, `hi`; `NULL` for `"dcdh"`); `meth`; `ctrl`; `reps`; and
#'   `level`.
#' @references
#'   Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-differences
#'   with multiple time periods. *Journal of Econometrics*, 225(2),
#'   200-230. \doi{10.1016/j.jeconom.2020.12.001}
#'
#'   Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
#'   difference-in-differences estimators. *Journal of Econometrics*,
#'   219(1), 101-122. \doi{10.1016/j.jeconom.2020.06.003}
#'
#'   Wooldridge, J. M. (2021). Two-way fixed effects, the two-way Mundlak
#'   regression, and difference-in-differences estimators.
#'   \doi{10.2139/ssrn.3906345}
#'
#'   de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed
#'   effects estimators with heterogeneous treatment effects. *American
#'   Economic Review*, 110(9), 2964-2996. \doi{10.1257/aer.20181169}
#' @examples
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
#' i <- match(d$town, LETTERS)
#' d$first <- c(2019, 2019, 2019, 2019, 2021, 2021, 2021, 2021,
#'              NA, NA, NA, NA)[i]
#' on <- !is.na(d$first) & d$year >= d$first
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * on + sin(seq_len(120))
#' s <- stgdid(d$rate, d$town, d$year, d$first, reps = 99, seed = 1)
#' s$overall
#' s$dynamic
#' # a by-law some towns repeal
#' d$law <- as.integer(on & !(i == 1 & d$year >= 2022))
#' stgdid(d$rate, d$town, d$year, trt = d$law, meth = "dcdh", reps = 49,
#'        seed = 1)$overall
#' @export
stgdid <- function(y, unit, time, first = NULL, trt = NULL, meth = "csa",
                   ctrl = "notyet", adj = NULL, est = "dr", clus = NULL,
                   reps = 499, seed = NULL, level = 0.95) {
  .Call(C_clab_r_stgdid, y, unit, time, first, trt, meth, ctrl, adj, est,
        clus, reps, seed, level)
}

#' @export
print.stgdid <- function(x, ...) {
  cat("stgdid: ", x[["meth"]],
      if (x[["meth"]] == "csa") paste0(", against units ", x[["ctrl"]]),
      ", ", x[["reps"]], " bootstrap samples\n", sep = "")
  print(x[["overall"]], row.names = FALSE)
  if (!is.null(x[["dynamic"]])) {
    cat("by event time:\n")
    print(x[["dynamic"]], row.names = FALSE)
  }
  invisible(x)
}
