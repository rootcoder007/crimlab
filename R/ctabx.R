#' Cross-tabulation with tests and measures of association
#'
#' Builds a two-way table from labels (counting rows, summing a count
#' column of an aggregate table, or counting distinct ids per cell) and
#' gives Pearson's chi-square test, the likelihood-ratio test, Cramer's V
#' with an interval and its bias-corrected form, the contingency
#' coefficient, phi for a 2 x 2 table, the share on the diagonal (the rows
#' whose label is also the column's, as a region at placement against the
#' region now), and on request Fisher's exact test with, for a 2 x 2 table,
#' the conditional odds ratio and its exact interval.
#'
#' Rows or columns with no counts are dropped before the statistics (their
#' number is `drop`), as an unused factor level carries no information;
#' the table `tab` keeps them. Cramer's V is `sqrt(chisq / (n (k - 1)))`,
#' `k` the smaller dimension, from the uncorrected statistic even when
#' `correct` is `TRUE`; its interval inverts the noncentral chi-square
#' distribution, `sqrt(ncp / (n (k - 1)))`, as the default of
#' `DescTools::CramerV()`; the bias-corrected V is Bergsma's. The exact
#' test enumerates every table with the observed margins and refuses one
#' that would take more than 1e7 tables. With `agg` other than `"n"`,
#' `"sum"` or `"dist"` the table holds values, not counts, and the
#' statistics are `NA`.
#'
#' @param row The row label of each record (a factor by its levels, text or
#'   numbers by their sorted values); or a whole table: a numeric matrix or
#'   2-d table of counts, with `col` `NULL`; or a `clpre` list from
#'   [clprep()].
#' @param col The column label of each record.
#' @param w `NULL`, or the value of each record that `agg` combines into
#'   its cell: a count (`"sum"`), a measure (`"max"`, `"min"`, `"mean"`), or
#'   an id (`"dist"`, the distinct ids per cell).
#' @param agg How the records of a cell combine: `"n"` (count them, the
#'   default without `w`), `"sum"` (the default with `w`), `"max"`, `"min"`,
#'   `"mean"` or `"dist"`. An empty cell is 0.
#' @param correct `TRUE` for the continuity correction of a 2 x 2
#'   chi-square, as [stats::chisq.test()]'s default; `FALSE` (here the
#'   default) for Pearson's statistic.
#' @param exact `TRUE` for Fisher's exact test.
#' @inheritParams cntreg
#' @return An object of class `ctabx`, a list: `tab`, the table (all labels);
#'   `expd` and `res`, the expected counts and Pearson residuals of the
#'   table without empty rows and columns; `stat`, a one-row data frame:
#'   `n`, `chisq`, `df`, `p`, `g2`, `g2p` (the likelihood-ratio statistic
#'   and its p-value), `v`, `vlo`, `vhi`, `vbc` (bias-corrected V), `cc`
#'   (the contingency coefficient), `phi`, `diag`, `fp` (Fisher's p-value),
#'   `or`, `orlo`, `orhi`; and `drop`, `agg`, `correct` and `level`.
#' @references
#'   Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley.
#'   ISBN 978-0-470-46363-5.
#'
#'   Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#'   *Journal of the Korean Statistical Society*, 42(3), 323-328.
#'   \doi{10.1016/j.jkss.2012.10.002}
#' @examples
#' # force type by weapon, one row per incident
#' wp <- c("none", "none", "knife", "gun", "none", "knife", "gun", "none")
#' fc <- c("hands", "hands", "CEW", "firearm", "CEW", "hands", "firearm",
#'         "hands")
#' ctabx(wp, fc)
#' # region at placement against region now: the diagonal share
#' ctabx(c("N", "N", "E", "W", "E"), c("N", "E", "E", "W", "E"))
#' # a published aggregate table: sum its count column
#' d <- data.frame(race = c("A", "A", "B", "B"), sex = c("F", "M", "F", "M"),
#'                 n = c(12, 40, 9, 51))
#' ctabx(d$race, d$sex, w = d$n, exact = TRUE)
#' @export
ctabx <- function(
  row,
  col = NULL,
  w = NULL,
  agg = NULL,
  correct = FALSE,
  exact = FALSE,
  level = 0.95
) {
  if (inherits(row, "clpre")) {
    return(.prep(ctabx, row, match.call(), parent.frame()))
  }
  .Call(C_clab_r_ctabx, row, col, w, agg, correct, exact, level)
}

#' @export
print.ctabx <- function(x, ...) {
  s <- x[["stat"]]
  cat(
    "ctabx: ",
    nrow(x[["tab"]]),
    " x ",
    ncol(x[["tab"]]),
    " table (",
    x[["agg"]],
    "), n = ",
    format(s[["n"]]),
    "\n",
    sep = ""
  )
  print(x[["tab"]])
  if (x[["drop"]] > 0L) {
    cat(x[["drop"]], "empty rows or columns left out of the statistics\n")
  }
  if (!is.na(s[["chisq"]])) {
    cat(
      "chi-square ",
      format(s[["chisq"]]),
      " on ",
      s[["df"]],
      " df, p = ",
      format(s[["p"]]),
      if (x[["correct"]]) " (continuity corrected)",
      "\nCramer's V ",
      format(s[["v"]]),
      " [",
      format(s[["vlo"]]),
      ", ",
      format(s[["vhi"]]),
      "], bias-corrected ",
      format(s[["vbc"]]),
      "\n",
      sep = ""
    )
  }
  if (!is.na(s[["diag"]])) {
    cat("diagonal share ", format(s[["diag"]]), "\n", sep = "")
  }
  if (!is.na(s[["fp"]])) {
    cat("Fisher's exact p = ", format(s[["fp"]]), "\n", sep = "")
  }
  invisible(x)
}
