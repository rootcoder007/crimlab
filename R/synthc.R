#' Synthetic control for one treated unit
#'
#' A weighted average of units never treated that tracks the treated unit
#' before the policy, as its counterfactual after; with the placebo test
#' that runs the same method on every comparison unit.
#'
#' The weights solve `min_w sum_{t < T0} (y_1t - sum_j w_j y_jt)^2` over
#' `w_j >= 0`, `sum_j w_j = 1`, the periods before adoption `T0`. The
#' problem is the nearest point of a convex hull to the origin, solved
#' exactly by Wolfe's minimum-norm-point algorithm. The effect in each
#' period is the gap `y_1t - sum_j w_j y_jt`. For inference each
#' comparison unit in turn is treated as if it had adopted, with the
#' others as its donors; the ratio of the root mean squared gap after
#' adoption to that before ranks the treated unit among the placebos, and
#' `p` is `(1 + the placebos with a ratio at least the treated unit's) /
#' (1 + the placebos)` (Abadie, Diamond and Hainmueller 2010). The weights
#' match the outcomes before adoption only, not separate predictors.
#'
#' @inheritParams evtst
#' @param time The period of each row: numbers; a balanced panel.
#' @param first The period the treated unit is first treated, on its rows;
#'   `NA`, 0 or `Inf` on the rows of the comparison units. Exactly one
#'   unit is treated, from a period after the first.
#' @param plac `TRUE` to run the placebo test.
#' @return An object of class `synthc`, a list: `weights` (`unit`, `w`);
#'   `path` (`time`, `treated`, `synth`, `gap`); `fit`, one row:
#'   `rmspepre`, `rmspepost`, their `ratio` and the placebo `p` (`NA`
#'   without the placebos); `plac` (`unit`, `ratio`; `NULL` without); and
#'   `unit`, the treated unit.
#' @references
#'   Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control
#'   methods for comparative case studies. *Journal of the American
#'   Statistical Association*, 105(490), 493-505.
#'   \doi{10.1198/jasa.2009.ap08746}
#'
#'   Wolfe, P. (1976). Finding the nearest point in a polytope.
#'   *Mathematical Programming*, 11, 128-149. \doi{10.1007/BF01580381}
#' @examples
#' d <- expand.grid(year = 2010:2024, town = LETTERS[1:8])
#' i <- match(d$town, LETTERS)
#' d$first <- ifelse(d$town == "A", 2019, NA)
#' d$rate <- 20 + i + (d$year - 2010) * (0.2 + i / 20) +
#'   sin(seq_len(120)) - 3 * (d$town == "A" & d$year >= 2019)
#' s <- synthc(d$rate, d$town, d$year, d$first)
#' s
#' s$path
#' @export
synthc <- function(y, unit, time, first, plac = TRUE) {
  .Call(C_clab_r_synthc, y, unit, time, first, plac)
}

#' @export
print.synthc <- function(x, ...) {
  f <- x[["fit"]]
  cat("synthc: unit ", x[["unit"]], ", root mean squared gap ",
      format(f[["rmspepre"]]), " before and ", format(f[["rmspepost"]]),
      " after", if (!is.na(f[["p"]])) paste0("; placebo p ",
                                             format(f[["p"]])),
      "\n", sep = "")
  w <- x[["weights"]]
  print(w[w[["w"]] > 0, ], row.names = FALSE)
  invisible(x)
}

#' Synthetic difference in differences
#'
#' The estimator of Arkhangelsky et al. (2021) for units adopting a policy
#' in one period: comparison units weighted so that their path before
#' adoption parallels the treated units', periods before adoption weighted
#' so that they predict the comparison units' level after, and the
#' difference in differences taken with both weights.
#'
#' With `N0` comparison units, `N1` treated, `T0` periods before adoption
#' and `T1` after, the unit weights solve `min (1/T0) sum_{t <= T0} (w0 +
#' sum_j w_j y_jt - ybar_1t)^2 + zeta^2 ||w||^2` over the simplex, with
#' `zeta = (N1 T1)^(1/4) sigma` and `sigma` the standard deviation of the
#' comparison units' first differences before adoption; the time weights
#' solve the same problem across comparison units, predicting each unit's
#' mean after adoption from its periods before, with `zeta = 1e-6 sigma`.
#' The estimate is `(ybar_1,post - sum_t l_t ybar_1t) - sum_j w_j
#' (ybar_j,post - sum_t l_t y_jt)`. The weights are the exact solutions of
#' these convex problems (Wolfe's minimum-norm point); the 'synthdid'
#' package runs Frank-Wolfe steps and then sparsifies, so its estimate
#' differs slightly. The standard error is the placebo one when one unit
#' is treated (a comparison unit drawn at random treated as if adopting,
#' the weights solved again, `reps` times; `sqrt((reps - 1) / reps)` times
#' the standard deviation), and otherwise the jackknife over units with
#' the weights held (Arkhangelsky et al. 2021, algorithms 3 and 4; `NA`
#' when one comparison unit carries all the weight).
#'
#' @inheritParams synthc
#' @param first The period each row's unit is first treated, one period for
#'   every treated unit; `NA`, 0 or `Inf` for the comparison units.
#' @param reps The placebo draws: a whole number, at least 2.
#' @param seed The seed of the placebo draws: `NULL` to draw it from R's
#'   generator, or a whole number from 0.
#' @inheritParams cbatt
#' @return An object of class `sdidx`, a list: `est` (`est`, `se`, `lo`,
#'   `hi`); `meth` (`"placebo"` or `"jackknife"`); `unitw` (`unit`, `w`);
#'   `timew` (`time`, `w`); `path` (`time`, `treated`, the treated mean,
#'   and `synth`, the weighted comparison units); `level`.
#' @references
#'   Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W. and
#'   Wager, S. (2021). Synthetic difference-in-differences. *American
#'   Economic Review*, 111(12), 4088-4118. \doi{10.1257/aer.20190159}
#' @examples
#' d <- expand.grid(year = 2010:2024, town = LETTERS[1:10])
#' i <- match(d$town, LETTERS)
#' d$first <- ifelse(i <= 3, 2019, NA)
#' d$rate <- 20 + i + (d$year - 2010) * (0.2 + i / 20) +
#'   sin(seq_len(150)) - 3 * (i <= 3 & d$year >= 2019)
#' sdidx(d$rate, d$town, d$year, d$first)
#' @export
sdidx <- function(y, unit, time, first, reps = 200, seed = NULL,
                  level = 0.95) {
  .Call(C_clab_r_sdidx, y, unit, time, first, reps, seed, level)
}

#' @export
print.sdidx <- function(x, ...) {
  cat("sdidx: standard error by ", x[["meth"]], "\n", sep = "")
  print(x[["est"]], row.names = FALSE)
  invisible(x)
}
