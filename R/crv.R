#' Cramer's V of a two-way table
#'
#' The strength of association in a contingency table as one number: 0 when
#' the rows and columns are independent in the sample, 1 when each row holds
#' a single column. It is the measure [ctabx()] reports beside its tests,
#' given here on its own for a table that is already counted.
#'
#' Cramer's V is `sqrt(chisq / (n (k - 1)))`, `k` the smaller dimension,
#' from Pearson's uncorrected statistic. In a small table it is biased
#' upward (it cannot be negative, so sampling noise alone gives a positive
#' value); with `biascor = TRUE`, the default, `v` is Bergsma's
#' bias-corrected form, and `vraw` holds the uncorrected one. The interval
#' `lo` to `hi` inverts the noncentral chi-square distribution and belongs
#' to the uncorrected V. Rows and columns with no counts are dropped
#' first, as an unused factor level carries no information.
#'
#' The p-value is that of the chi-square approximation, which rests on the
#' expected counts: `minexp` is the smallest of them, and below about 5
#' the exact test of [ctabx()] (`exact = TRUE`) is the one to report.
#'
#' Errors: `tbl` not a numeric matrix or two-way table of finite,
#' non-negative counts; `biascor` not `TRUE` or `FALSE`; a `level` outside
#' (0, 1).
#'
#' @param tbl A two-way table of counts: a numeric matrix or a 2-d table.
#' @param biascor `TRUE` (the default) for Bergsma's bias-corrected V in
#'   `v`; `FALSE` for the uncorrected one.
#' @inheritParams ctabx
#' @return A named numeric vector of class `cramerv`: `v`, `vraw`, `lo`,
#'   `hi`, `chisq`, `df`, `p`, `n`, `minexp` (the smallest expected count)
#'   and `drop` (the empty rows and columns left out). In a table with
#'   fewer than two rows or columns that hold counts the statistics are
#'   `NA`.
#' @references
#'   Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton
#'   University Press.
#'
#'   Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#'   *Journal of the Korean Statistical Society*, 42(3), 323-328.
#'   \doi{10.1016/j.jkss.2012.10.002}
#' @seealso [ctabx()] for the table from records, the tests and the exact
#'   test.
#' @examples
#' # use of force by weapon, counted
#' tab <- matrix(c(30, 6, 2, 9, 14, 3, 1, 4, 11), 3,
#'               dimnames = list(c("none", "knife", "gun"),
#'                               c("hands", "CEW", "firearm")))
#' cramerv(tab)
#' cramerv(tab, biascor = FALSE)[["v"]]
#' # independence: the corrected V is 0
#' cramerv(matrix(c(10, 20, 10, 20), 2))[["v"]]
#' try(cramerv(1:4))
#' @export
cramerv <- function(tbl, biascor = TRUE, level = 0.95) {
  fn <- "cramerv"
  biascor <- .ckflag(biascor, "biascor", fn)
  if (
    !(is.matrix(tbl) || inherits(tbl, "table")) ||
      !is.numeric(tbl) ||
      length(dim(tbl)) != 2L
  ) {
    .ckstop(
      fn,
      paste(
        "`tbl` must be a two-way table of counts",
        "(a numeric matrix or a 2-d table)"
      )
    )
  }
  x <- .clas(fn, ctabx(tbl, level = level))
  s <- x[["stat"]]
  e <- x[["expd"]]
  out <- c(
    v = if (biascor) s[["vbc"]] else s[["v"]],
    vraw = s[["v"]],
    lo = s[["vlo"]],
    hi = s[["vhi"]],
    chisq = s[["chisq"]],
    df = s[["df"]],
    p = s[["p"]],
    n = s[["n"]],
    minexp = min(e),
    drop = x[["drop"]]
  )
  structure(
    out,
    biascor = biascor,
    level = level,
    class = c("cramerv", "numeric")
  )
}

#' @export
print.cramerv <- function(x, ...) {
  v <- unclass(x)
  cat(
    "cramerv: V = ",
    format(v[["v"]]),
    if (attr(x, "biascor")) {
      " (bias-corrected; uncorrected "
    } else {
      " (uncorrected; "
    },
    if (attr(x, "biascor")) format(v[["vraw"]]),
    if (attr(x, "biascor")) ", " else "",
    format(100 * attr(x, "level")),
    "% interval ",
    format(v[["lo"]]),
    " to ",
    format(v[["hi"]]),
    ")\nchi-square ",
    format(v[["chisq"]]),
    " on ",
    format(v[["df"]]),
    " df, p = ",
    format(v[["p"]]),
    ", n = ",
    format(v[["n"]]),
    "\n",
    sep = ""
  )
  if (!is.na(v[["minexp"]]) && v[["minexp"]] < 5) {
    cat(
      "smallest expected count ",
      format(v[["minexp"]]),
      ": report the exact test, ctabx(exact = TRUE)\n",
      sep = ""
    )
  }
  if (v[["drop"]] > 0) {
    cat(format(v[["drop"]]), "empty rows or columns left out\n")
  }
  invisible(x)
}
