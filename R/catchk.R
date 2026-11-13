#' Audit categorical inputs before an analysis
#'
#' Reports, for each variable, what the functions of this package would
#' make of it: its storage type, its levels and their counts, missing
#' values, unused levels, levels counted fewer than `nmin` times, labels
#' that differ only by case or white space, a 0/1 or logical coding,
#' integer codes that look like categories, value labels (as haven makes)
#' that do not match the values, and which level a default would take as
#' the reference, and why. Nothing is changed.
#'
#' Levels are those of a factor, in level order; `FALSE` and `TRUE` for a
#' logical, and `0` and `1` for numbers that are all 0 or 1, both always
#' listed; otherwise the distinct values as text, in order of first
#' appearance, and for a labelled vector then the labelled values that do
#' not occur. Doubles that are not all whole numbers, dates, times, lists
#' and matrices are not categorical: they are reported with their type and
#' missing values only.
#'
#' The default reference: `FALSE` for a logical and `0` for 0/1 numbers
#' (`sntbnd()` takes `TRUE` and 1 as the treatment), the first level of a
#' factor (`ref = 1` on a table of it), and none otherwise: character
#' values and codes have no order the package uses, so the reference must
#' be stated, and a labelled vector is refused until it is decoded.
#' Labels are compared after trimming white space, making each run of it
#' one space and ASCII letters lower case; other letters are compared as
#' they are.
#'
#' @param x A vector, a data frame, or a list of vectors (each one a
#'   variable, named by its name or position).
#' @param ref `NULL`, or the reference level the analysis intends: one
#'   level for every variable, or levels named by variable (a named vector
#'   or list). Matched exactly to the levels as text; a factor by its
#'   label.
#' @param nmin Levels counted fewer times (but at least once) are reported
#'   as low: one whole number from 1.
#' @return An object of class `catchk`, a list: `ok` (`TRUE` when no
#'   variable has a problem), `nmin`, `vars`, a data frame with one row per
#'   variable, and `levels`, one row per level. `vars`: `var`, `type`, `n`
#'   (values), `nna` (missing), `nlev` (levels), `unused` (levels never
#'   seen), `low` (levels seen fewer than `nmin` times), `twin` (levels
#'   with a label that differs from another only by case or white space),
#'   `bin` (logical or 0/1), `codes` (whole numbers, not 0/1, with at most
#'   10 distinct values: likely category codes), `nolab` (values seen with
#'   no value label), `dlab` (value labels shared by two values), `ref` (the
#'   default reference, `NA` for none), `why`, `want` (the stated
#'   reference), `state` (`"default"`: the stated level is the default;
#'   `"other"`: present, but must be passed explicitly; `"low"`, `"empty"`
#'   (never seen), `"case"` (present only up to case or white space),
#'   `"absent"`; `NA` when none is stated) and `ok`. `nlev`, `unused`,
#'   `low` and `twin` are `NA` for a variable that is not categorical.
#'   `levels`: `var`, `level`, `label` (the value label, `NA` where none),
#'   `n`, `low`, `twin` (the level it differs from only by case or white
#'   space, `NA` where none) and `ref` (`TRUE` for the default reference).
#'   A variable has a problem when it has missing values, unused or low
#'   levels, twins, value labels, likely codes, fewer than two levels seen,
#'   or a stated reference in state `"low"`, `"empty"`, `"case"` or
#'   `"absent"`.
#' @examples
#' d <- data.frame(sex = c("Male", "male", "Female", "Female", NA),
#'                 custody = c(TRUE, FALSE, TRUE, TRUE, FALSE),
#'                 court = c(1L, 2L, 2L, 3L, 1L))
#' catchk(d, nmin = 2)
#' catchk(factor(c("b", "a", "b"), levels = c("a", "b", "c")), ref = "b")
#' @export
catchk <- function(x, ref = NULL, nmin = 5) {
  .Call(C_clab_r_catchk, x, ref, nmin)
}

#' @export
print.catchk <- function(x, ...) {
  v <- x[["vars"]]
  l <- x[["levels"]]
  q <- function(s) paste0("\"", s, "\"", collapse = ", ")
  bad <- c(low = "seen fewer than nmin times", empty = "never seen",
           case = "present only with other case or spacing",
           absent = "not a level", default = "", other = "")
  p <- vapply(seq_len(nrow(v)), function(i) {
    w <- as.list(v[i, ])
    li <- l[l[["var"]] == w[["var"]], ]
    # the levels of one fold together, keyed by the first of them
    tw <- li[!is.na(li[["twin"]]), ]
    at <- match(tw[["twin"]], li[["level"]]) > match(tw[["level"]],
                                                     li[["level"]])
    root <- ifelse(at, tw[["level"]], tw[["twin"]])
    tw <- vapply(split(tw[["level"]], factor(root, unique(root))), q, "")
    # no stated reference (NA) reads as ""
    st <- c(bad, "")[[match(w[["state"]], names(bad), length(bad) + 1L)]]
    s <- c(paste(w[["nna"]], "missing"),
           paste("never seen:", q(li[["level"]][li[["n"]] == 0L])),
           paste0("seen fewer than ", x[["nmin"]], " times: ",
                  q(li[["level"]][li[["low"]]])),
           paste("differ only by case or spaces:", paste(tw, collapse = "; ")),
           "value labels: decode it first",
           paste(w[["nolab"]], "values with no label"),
           paste(w[["dlab"]], "labels shared by two values"),
           "numbers that look like category codes",
           "fewer than two levels seen",
           paste0("reference ", q(w[["want"]]), ": ", st))
    k <- c(w[["nna"]] > 0L, isTRUE(w[["unused"]] > 0L),
           isTRUE(w[["low"]] > 0L), length(tw) > 0L,
           w[["type"]] == "labelled", w[["nolab"]] > 0L, w[["dlab"]] > 0L,
           w[["codes"]], isTRUE(w[["nlev"]] - w[["unused"]] < 2L),
           nzchar(st))
    paste(s[k], collapse = "; ")
  }, "")
  cat("catchk: ", nrow(v), if (nrow(v) == 1L) " variable" else " variables",
      ", ", sum(!v[["ok"]]), " with problems\n", sep = "")
  print(data.frame(var = v[["var"]], type = v[["type"]], nlev = v[["nlev"]],
                   ref = v[["ref"]], why = v[["why"]]), row.names = FALSE)
  for (i in which(nzchar(p))) cat(v[["var"]][i], ": ", p[i], "\n", sep = "")
  invisible(x)
}
