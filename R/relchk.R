#' Which labelling reproduces a reported contrast?
#'
#' A reported odds ratio, risk ratio or risk difference depends on which
#' group is the reference, which outcome is the event and which margin of
#' the table holds the groups. `relchk()` recomputes the contrast from the
#' table under every such reading and says which ones reproduce the
#' reported value: the labelling the report states, another one, or none.
#'
#' A reading compares group `group` with group `ref` for outcome `event`,
#' the other outcomes pooled as non-events: with `a` and `b` the events and
#' non-events of `group`, `c` and `d` those of `ref`, the odds ratio is
#' `a d / (b c)`, the risk ratio `(a / (a + b)) / (c / (c + d))` and the
#' risk difference `a / (a + b) - c / (c + d)`. Every ordered pair of
#' groups and every outcome is read, with the groups in the rows and, the
#' table transposed, in the columns, and all three kinds of contrast.
#'
#' No correction is made for zero cells: a reading whose denominator is 0
#' is `Inf` (a positive numerator) or `NaN` (0 over 0) and matches nothing.
#' Counts need not be whole numbers (weighted totals are read as they are).
#'
#' For a 2 x 2 table the odds ratio inverts when the reference group or the
#' event changes and is unchanged when both change or the table is
#' transposed; the risk ratio inverts with the reference group, and the
#' risk difference changes sign with either; `proof("relchk")` names the
#' Lean theorems.
#'
#' @param x The table of counts, groups in rows and outcomes in columns: a
#'   matrix, a 2-d table, a data frame of numeric columns, or a list of
#'   numeric vectors of one length (the columns). At least 2 x 2; counts
#'   finite and at least 0 with no NA, no row or column all 0, and no label
#'   missing or given twice. Rows and columns without names are labelled
#'   "1", "2", ... .
#' @param value The reported value: one finite number.
#' @param type The kind of contrast the report claims: `"or"` (odds ratio),
#'   `"rr"` (risk ratio) or `"rd"` (risk difference).
#' @param group,ref,event The labelling the report claims, each `NULL` (not
#'   stated) or a label (matched exactly, case included; a factor by its
#'   label) or a position: `group` and `ref` name rows of `x`, `event` a
#'   column.
#' @param digits The decimals `value` was rounded to: a reading matches
#'   when it lies within half a unit of the last decimal, at either end, so
#'   that both rounding conventions at the half are accepted. By default
#'   the fewest decimals that print `value` exactly (give `digits = 2` for
#'   a value reported as 2.30). Give `digits` or `tol`, not both.
#' @param tol Instead of `digits`, the largest difference allowed: one
#'   number from 0. Either way, a reading also matches when it differs only
#'   in the last few bits.
#' @return An object of class `relchk`, a list: `verdict` (`"holds"`: a
#'   reading the claim allows matches; `"other"`: only readings it does not
#'   allow match; `"none"`), `ok` (`TRUE` when the claim holds), `value`,
#'   `tol` (the difference allowed), `digits` (`NA` when `tol` was given or
#'   no decimals print `value` exactly; then `tol` is 0), `type`, `claim`
#'   (the stated `group`, `ref` and `event`, `NA` where not stated) and
#'   `reads`, a data frame with one row per reading: `grp` (`"rows"` or
#'   `"cols"`: the margin holding the groups), `type`, `group`, `ref`,
#'   `event`, `value`, `match` and `claim` (`TRUE` for the readings the
#'   claim allows: groups in rows, the kind `type`, and the stated labels).
#' @examples
#' tab <- matrix(c(20, 10, 80, 90), 2,
#'               dimnames = list(c("A", "B"), c("yes", "no")))
#' # 2.25 is the odds ratio of A against B for "yes"
#' relchk(tab, 2.25, ref = "B", event = "yes")
#' # 0.44 is the same contrast read with the reference group swapped
#' relchk(tab, 0.44, ref = "B", event = "yes")
#' # a risk ratio reported as an odds ratio
#' relchk(tab, 2.0, type = "or", digits = 1)
#' @export
relchk <- function(x, value, type = "or", group = NULL, ref = NULL,
                   event = NULL, digits = NULL, tol = NULL) {
  .Call(C_clab_r_relchk, x, value, type, group, ref, event, digits, tol)
}

#' @srrstats {EA4.1} The print methods show reported values with
#'   format(), at full precision unless options(digits) says otherwise;
#'   the numbers in the result are unrounded.
#' @srrstats {EA5.2} print.relchk() and print.catchk() print the verdict
#'   and the readings or problems that explain it; the result holds all
#'   of it.
#' @export
print.relchk <- function(x, ...) {
  r <- x[["reads"]]
  cl <- x[["claim"]]
  st <- !is.na(cl)
  d <- x[["digits"]]
  cat("relchk: reported ", format(x[["value"]]), ", ",
      if (is.na(d)) {
        paste("within", format(x[["tol"]]))
      } else {
        paste("to", d, if (d == 1L) "decimal" else "decimals")
      }, "\n", sep = "")
  cat("claim: ", x[["type"]],
      if (any(st)) paste0(", ", names(cl)[st], " ", cl[st], collapse = ""),
      "\n", sep = "")
  cat("verdict: ",
      switch(x[["verdict"]],
             holds = "the claim holds",
             other = "only readings the claim does not allow match",
             none = "no reading matches"), "\n", sep = "")
  m <- r[r[["match"]], c("grp", "type", "group", "ref", "event", "value",
                         "claim")]
  if (nrow(m) > 0L) print(m, row.names = FALSE)
  nu <- sum(!is.finite(r[["value"]]))
  if (nu > 0L) {
    cat(nu, " of ", nrow(r), " readings undefined: a zero denominator\n",
        sep = "")
  }
  invisible(x)
}
