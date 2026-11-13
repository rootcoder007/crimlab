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
