# Rules a table must satisfy, declared by the analyst and run before an
# analysis. catchk() audits categories and relchk() a reported contrast;
# these say what only the analyst knows: that an age is not negative, that
# an identifier is unique, that a count never falls.
#
# A ready-made rule is data: its kind and what it was given, with no
# function in it, so that the same rule runs in every language of the
# package (cli/jval.inc holds the C++ of valrules()). Only a rule made by
# rul(), which holds an R function, is R's alone.

#' Declare rules a table must satisfy, and check them
#'
#' `rul()` declares one rule: a name, a function that says whether the
#' data pass, the column it reads (or the whole table), how serious a
#' failure is and the message to give. `valrules()` runs a list of rules
#' on a data frame and returns the ones that failed, with the rows that
#' broke them. The ready-made rules ([ruleinset()] and the others) cover
#' the usual checks, so that they are written once and not by hand.
#'
#' The function `expr` receives the column `col` as a vector, or the whole
#' data frame when `col` is `NULL`. It returns either one `TRUE` or
#' `FALSE` for the table, or one value per row; in the row-wise form a
#' missing answer counts as a failure, since a rule that cannot decide has
#' not been met. A rule for a column the table does not have is skipped. A
#' rule whose function stops with an error, or returns something that is
#' not logical, is reported as a failure that says so: a broken rule never
#' passes silently. Two rules of one name leave one entry, the later one.
#'
#' A ready-made rule holds no function: it is a list of its kind (`rule`)
#' and what it was given, and `valrules()` also takes it written by hand,
#' as `list(rule = "btw", col = "year", lo = 2000, hi = 2026)`. That form
#' is the one the other languages of the package use, where `valrules()`
#' and the ready-made rules exist too; a rule made by `rul()` holds an R
#' function and runs in R only.
#'
#' Errors: a `name` or `col` that is not one non-empty string; `expr` not
#' a function; `sev` not `"warning"` or `"fatal"`; `data` not a data
#' frame; `rules` not a rule or a list of rules.
#'
#' @param name The rule's name: one non-empty string. It names the rule in
#'   the result.
#' @param expr A function of one argument that returns `TRUE`/`FALSE` for
#'   the table or one logical value per row.
#' @param col The column the rule reads, or `NULL` (the default) for a
#'   rule on the whole data frame.
#' @param sev `"warning"` (the default) or `"fatal"`: how serious a
#'   failure is. `valrules()` records it; the caller decides what follows.
#' @param msg The message of a failure, or `NULL` for "Rule '<name>'
#'   failed".
#' @param data The table the rules are run on: a data frame.
#' @param rules One rule or a list of rules.
#' @return `rul()` an object of class `clrule`, a list: `rule` (`"fn"`),
#'   `name`, `col` (when given), `sev`, `msg` and `expr`. `valrules()` an
#'   object of class `valrules`: a list with one element per failed rule,
#'   named by the rule, each a list of `severity`, `message` and, for a
#'   row-wise rule, `nfail` (the number of failing rows) and `rows` (the
#'   first ten of them). A table that passes every rule gives an empty
#'   list.
#' @references
#'   van der Loo, M. P. J. and de Jonge, E. (2021). Data validation
#'   infrastructure for R. *Journal of Statistical Software*, 97(10), 1-31.
#'   \doi{10.18637/jss.v097.i10}
#'
#'   de Jonge, E. and van der Loo, M. (2018). *Statistical Data Cleaning
#'   with Applications in R*. Wiley. \doi{10.1002/9781118897126}
#' @seealso [ruleinset()] for the ready-made rules; [catchk()] and
#'   [prechk()] for other checks before an analysis.
#' @examples
#' d <- data.frame(id = c(1, 2, 2, 4), age = c(34, -1, 51, NA),
#'                 charge = c("assault", "theft", "theft", "arson"))
#' r <- list(
#'   rul("age_not_negative", function(v) is.na(v) | v >= 0, col = "age",
#'       msg = "an age is negative"),
#'   ruleuniq("id", sev = "fatal"),
#'   rulenull("age"),
#'   ruleinset("charge", c("assault", "theft", "drugs"))
#' )
#' v <- valrules(d, r)
#' v
#' v[["id_unique"]][["rows"]]
#' # a table that passes
#' length(valrules(d[1, ], r))
#' # a rule written as data
#' valrules(d, list(rule = "btw", col = "age", lo = 0, hi = 120))
#' try(rul("", function(v) TRUE))
#' @export
rul <- function(
  name,
  expr,
  col = NULL,
  sev = c("warning", "fatal"),
  msg = NULL
) {
  fn <- "rul"
  sev <- .ckmarg("sev", fn)
  name <- .rlstr(name, "name", fn)
  if (!is.function(expr)) {
    .ckstop(fn, "`expr` must be a function")
  }
  if (!is.null(col)) {
    col <- .rlstr(col, "col", fn)
  }
  msg <- if (is.null(msg)) {
    sprintf("Rule '%s' failed", name)
  } else {
    .rlstr(msg, "msg", fn)
  }
  .rlnew("fn", name, col, sev, msg, list(expr = expr))
}

#' @rdname rul
#' @export
valrules <- function(data, rules) {
  fn <- "valrules"
  if (!is.data.frame(data)) {
    .ckstop(fn, "`data` must be a data frame")
  }
  one <- function(r) is.list(r) && is.character(r[["rule"]])
  if (one(rules)) {
    rules <- list(rules)
  }
  if (!is.list(rules) || !all(vapply(rules, one, NA))) {
    .ckstop(
      fn,
      paste(
        "`rules` must be a rule or a list of rules,",
        "from rul() or the rule functions"
      )
    )
  }
  out <- list()
  for (r in rules) {
    r <- .rlspec(r, fn)
    iss <- .rlrun(r, data)
    if (!is.null(iss)) {
      out[[r[["name"]]]] <- iss
    }
  }
  structure(out, nrule = length(rules), class = "valrules")
}

#' @export
print.clrule <- function(x, ...) {
  cat(
    "rule ",
    x[["name"]],
    " (",
    x[["sev"]],
    ", ",
    if (is.null(x[["col"]])) {
      "on the table"
    } else {
      paste0("on `", x[["col"]], "`")
    },
    "): ",
    x[["msg"]],
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
print.valrules <- function(x, ...) {
  n <- length(x)
  cat(
    "valrules: ",
    attr(x, "nrule"),
    " rules, ",
    n,
    " failed",
    if (n) {
      paste0(
        " (",
        sum(vapply(x, function(i) i[["severity"]] == "fatal", NA)),
        " fatal)"
      )
    },
    "\n",
    sep = ""
  )
  for (k in names(x)) {
    i <- x[[k]]
    cat(
      "  ",
      i[["severity"]],
      "  ",
      k,
      ": ",
      i[["message"]],
      if (!is.null(i[["rows"]])) {
        paste0(
          "; rows ",
          toString(i[["rows"]]),
          if (i[["nfail"]] > length(i[["rows"]])) ", ..."
        )
      },
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Ready-made rules
#'
#' The checks that account for most of what a table is asked: a value in
#' a set, a number in a range, no missing value, no repeated identifier, a
#' pattern, an order, an outlier, complete rows, distinct rows, a number
#' of columns. Each returns a rule for [valrules()], as [rul()] does.
#'
#' `ruleinset()`: every value is one of `set` (compared as text).
#' `rulebtw()`: every value lies in `[lo, hi]`. `rulenull()`: no value is
#' missing. `ruleuniq()`: no value repeats; the rows reported are the
#' repeats, not the first occurrence. `rulerx()`: every value matches the
#' regular expression `pat`. `ruleincr()`: the column never falls (or,
#' with `strict`, always rises) down the rows; missing values are passed
#' over. `rulemad()`: every value lies within `n` median absolute
#' deviations of the median, the deviation scaled by 1.4826 as
#' [stats::mad()] scales it (so that it estimates the standard deviation
#' of normal data; Leys et al. 2013 argue for this screen over one built
#' on the mean and the standard deviation, which the outliers themselves
#' inflate); when the deviation is 0 (more than half the values equal)
#' any other value fails. `rulecmpl()`: no row has a missing value.
#' `ruledstn()`: no row repeats, over all columns or over `cols`.
#' `rulecolcnt()`: the table has `n` columns.
#'
#' With `napass = TRUE`, the default, a missing value passes a rule about
#' values (missingness is `rulenull()`'s question); with `FALSE` it fails.
#' `rulebtw()`, `ruleincr()` and `rulemad()` read numbers (or dates, or
#' `TRUE` and `FALSE`): on a column of text they are reported as rules
#' that could not be evaluated, not passed.
#'
#' A pattern is read by [grepl()] with `perl = TRUE`; the other languages
#' of the package read it as ECMAScript, and the syntax the two share
#' (classes, anchors, counts, alternation, groups) means the same in
#' both. A pattern is at most 256 bytes, and a value it is tried on at
#' most 4096.
#'
#' Errors: `col` not one non-empty string; `lo` or `hi` not one number, or
#' `lo` above `hi`; `n` not one positive number (`rulecolcnt()`: one whole
#' number of at least 0); `pat` not one string of at most 256 bytes, or
#' not a regular expression; `cols` neither `NULL` nor column names; `sev`
#' not `"warning"` or `"fatal"`.
#'
#' @param col The column the rule reads: one non-empty string.
#' @param set The allowed values.
#' @param lo,hi The bounds, both included.
#' @param pat A regular expression.
#' @param strict `TRUE` for strictly increasing.
#' @param n `rulemad()`: the number of median absolute deviations allowed.
#'   `rulecolcnt()`: the number of columns.
#' @param cols The columns over which rows must be distinct, or `NULL`
#'   (all of them).
#' @param napass `TRUE` (the default) lets a missing value pass.
#' @inheritParams rul
#' @return A rule: an object of class `clrule`, a list of `rule` (the
#'   kind), `name`, `col` (for a rule on one column), `sev`, `msg` and
#'   what the rule was given.
#' @references
#'   Leys, C., Ley, C., Klein, O., Bernard, P. and Licata, L. (2013).
#'   Detecting outliers: do not use standard deviation around the mean,
#'   use absolute deviation around the median. *Journal of Experimental
#'   Social Psychology*, 49(4), 764-766. \doi{10.1016/j.jesp.2013.03.013}
#'
#'   Hampel, F. R. (1974). The influence curve and its role in robust
#'   estimation. *Journal of the American Statistical Association*,
#'   69(346), 383-393. \doi{10.1080/01621459.1974.10482962}
#' @seealso [rul()] and [valrules()].
#' @examples
#' d <- data.frame(year = c(2019, 2020, 2021, 2021),
#'                 n = c(410, 395, 9000, 402),
#'                 div = c("11", "14", "14", "5X"))
#' valrules(d, list(
#'   rulebtw("year", 2000, 2026),
#'   ruleincr("year", strict = TRUE),
#'   rulemad("n", 5),
#'   rulerx("div", "^[0-9]{2}$"),
#'   ruledstn(c("year", "div")),
#'   rulecolcnt(3),
#'   rulecmpl()
#' ))
#' try(rulebtw("year", 2026, 2000))
#' @export
ruleinset <- function(
  col,
  set,
  napass = TRUE,
  sev = c("warning", "fatal")
) {
  fn <- "ruleinset"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  napass <- .ckflag(napass, "napass", fn)
  if (!is.atomic(set) || !length(set)) {
    .ckstop(fn, "`set` must be a vector of at least one allowed value")
  }
  .rlnew(
    "inset",
    paste0(col, "_in_set"),
    col,
    sev,
    sprintf("Column '%s' has values outside the allowed set", col),
    list(set = as.character(set), napass = napass)
  )
}

#' @rdname ruleinset
#' @export
rulebtw <- function(
  col,
  lo,
  hi,
  napass = TRUE,
  sev = c("warning", "fatal")
) {
  fn <- "rulebtw"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  napass <- .ckflag(napass, "napass", fn)
  lo <- .rlone(lo, "lo", fn)
  hi <- .rlone(hi, "hi", fn)
  if (lo > hi) {
    .ckstop(fn, "`lo` must not be above `hi`")
  }
  .rlnew(
    "btw",
    paste0(col, "_between"),
    col,
    sev,
    sprintf(
      "Column '%s' has values outside [%s, %s]",
      col,
      as.character(lo),
      as.character(hi)
    ),
    list(lo = lo, hi = hi, napass = napass)
  )
}

#' @rdname ruleinset
#' @export
rulenull <- function(col, sev = c("warning", "fatal")) {
  fn <- "rulenull"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  .rlnew(
    "null",
    paste0(col, "_not_null"),
    col,
    sev,
    sprintf("Column '%s' has missing values", col)
  )
}

#' @rdname ruleinset
#' @export
ruleuniq <- function(col, sev = c("warning", "fatal")) {
  fn <- "ruleuniq"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  .rlnew(
    "uniq",
    paste0(col, "_unique"),
    col,
    sev,
    sprintf("Column '%s' has duplicate values", col)
  )
}

#' @rdname ruleinset
#' @export
rulerx <- function(col, pat, napass = TRUE, sev = c("warning", "fatal")) {
  fn <- "rulerx"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  napass <- .ckflag(napass, "napass", fn)
  pat <- .rlstr(pat, "pat", fn)
  if (nchar(pat, "bytes") > 256L) {
    .ckstop(fn, "`pat` must be at most 256 bytes")
  }
  ok <- tryCatch(
    is.logical(suppressWarnings(grepl(pat, "", perl = TRUE))),
    error = function(e) FALSE
  )
  if (!ok) {
    .ckstop(fn, "`pat` is not a valid regular expression")
  }
  .rlnew(
    "rx",
    paste0(col, "_regex"),
    col,
    sev,
    sprintf("Column '%s' has values not matching /%s/", col, pat),
    list(pat = pat, napass = napass)
  )
}

#' @rdname ruleinset
#' @export
ruleincr <- function(col, strict = FALSE, sev = c("warning", "fatal")) {
  fn <- "ruleincr"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  strict <- .ckflag(strict, "strict", fn)
  .rlnew(
    "incr",
    paste0(col, "_increasing"),
    col,
    sev,
    sprintf(
      "Column '%s' is not %s",
      col,
      if (strict) "strictly increasing" else "non-decreasing"
    ),
    list(strict = strict)
  )
}

#' @rdname ruleinset
#' @export
rulemad <- function(col, n = 3, napass = TRUE, sev = c("warning", "fatal")) {
  fn <- "rulemad"
  col <- .rlstr(col, "col", fn)
  sev <- .ckmarg("sev", fn)
  napass <- .ckflag(napass, "napass", fn)
  n <- .rlone(n, "n", fn)
  if (n <= 0) {
    .ckstop(fn, "`n` must be positive")
  }
  .rlnew(
    "mad",
    paste0(col, "_within_mads"),
    col,
    sev,
    sprintf(
      "Column '%s' has values more than %s MADs from its median",
      col,
      as.character(n)
    ),
    list(n = n, napass = napass)
  )
}

#' @rdname ruleinset
#' @export
rulecmpl <- function(sev = c("warning", "fatal")) {
  sev <- .ckmarg("sev", "rulecmpl")
  .rlnew("cmpl", "complete_rows", NULL, sev, "Some rows have missing values")
}

#' @rdname ruleinset
#' @export
ruledstn <- function(cols = NULL, sev = c("warning", "fatal")) {
  fn <- "ruledstn"
  sev <- .ckmarg("sev", fn)
  if (
    !is.null(cols) &&
      !(is.character(cols) && length(cols) && !anyNA(cols) && all(nzchar(cols)))
  ) {
    .ckstop(fn, "`cols` must be NULL or a character vector of column names")
  }
  .rlnew(
    "dstn",
    "distinct_rows",
    NULL,
    sev,
    if (is.null(cols)) {
      "Some rows are duplicated"
    } else {
      sprintf("Some rows are duplicated on (%s)", toString(cols))
    },
    if (!is.null(cols)) list(cols = cols)
  )
}

#' @rdname ruleinset
#' @export
rulecolcnt <- function(n, sev = c("warning", "fatal")) {
  fn <- "rulecolcnt"
  sev <- .ckmarg("sev", fn)
  n <- .rlone(n, "n", fn)
  if (n < 0 || n != floor(n)) {
    .ckstop(fn, "`n` must be a whole number of at least 0")
  }
  .rlnew(
    "colcnt",
    "col_count",
    NULL,
    sev,
    sprintf("Expected %.0f columns", n),
    list(n = n)
  )
}

#' A rule: its kind, name, column (when it reads one), severity and
#' message, then what it was given.
#' @noRd
.rlnew <- function(rule, name, col, sev, msg, more = NULL) {
  structure(
    c(
      list(rule = rule, name = name),
      if (!is.null(col)) list(col = col),
      list(sev = sev, msg = msg),
      more
    ),
    class = "clrule"
  )
}

#' A rule as valrules() runs it: a rule of this package as it is, a rule
#' written as data made by the function of its kind (which checks it).
#' @noRd
.rlspec <- function(r, fn) {
  if (inherits(r, "clrule")) {
    return(r)
  }
  make <- list(
    inset = ruleinset,
    btw = rulebtw,
    null = rulenull,
    uniq = ruleuniq,
    rx = rulerx,
    incr = ruleincr,
    mad = rulemad,
    cmpl = rulecmpl,
    dstn = ruledstn,
    colcnt = rulecolcnt
  )
  kind <- r[["rule"]]
  if (length(kind) != 1L || is.na(kind) || !kind %in% names(make)) {
    .ckstop(
      fn,
      paste(
        "`rules` must be a rule or a list of rules,",
        "from rul() or the rule functions"
      )
    )
  }
  f <- make[[kind]]
  a <- r[intersect(names(r), names(formals(f)))]
  # what the rule needs and was not given is NULL, so that the function
  # refuses it in its own words
  for (k in setdiff(
    names(formals(f)),
    c(names(a), "napass", "sev", "strict", "cols")
  )) {
    a[k] <- list(NULL)
  }
  if (kind == "mad" && is.null(a[["n"]])) {
    a[["n"]] <- NULL
  }
  do.call(f, a)
}

#' One non-empty string that is valid text, or the refusal of fn.
#' @noRd
.rlstr <- function(x, what, fn) {
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x) ||
      !validEnc(x)
  ) {
    .ckstop(fn, "`%s` must be one non-empty string", what)
  }
  x
}

#' One finite number, or the refusal of fn.
#' @noRd
.rlone <- function(x, what, fn) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    .ckstop(fn, "`%s` must be one number", what)
  }
  x
}

#' The row-wise answers ok with the missing values of v passed (napass) or
#' failed.
#' @noRd
.rlna <- function(ok, v, napass) {
  if (napass) ok | is.na(v) else ok & !is.na(v)
}

#' The numbers of a column a rule about numbers reads, or an error.
#' @noRd
.rlnum <- function(v) {
  if (!(is.numeric(v) || is.logical(v) || inherits(v, c("Date", "POSIXt")))) {
    stop("the column is not numeric", call. = FALSE)
  }
  as.numeric(v)
}

#' The answer of a ready-made rule r on its column v (or on the table d).
#' @noRd
.rlans <- function(r, v) {
  napass <- r[["napass"]]
  switch(
    r[["rule"]],
    fn = r[["expr"]](v),
    inset = .rlna(as.character(v) %in% r[["set"]], v, napass),
    btw = {
      x <- .rlnum(v)
      .rlna(x >= r[["lo"]] & x <= r[["hi"]], x, napass)
    },
    null = !is.na(v),
    uniq = !duplicated(v),
    rx = {
      s <- as.character(v)
      if (any(nchar(s, "bytes") > 4096L, na.rm = TRUE)) {
        stop("a value is longer than 4096 bytes", call. = FALSE)
      }
      .rlna(grepl(r[["pat"]], s, perl = TRUE), v, napass)
    },
    incr = {
      x <- .rlnum(v)
      d <- diff(x[!is.na(x)])
      if (r[["strict"]]) all(d > 0) else all(d >= 0)
    },
    mad = {
      x <- .rlnum(v)
      k <- !is.na(x)
      if (!any(k)) {
        return(rep(napass, length(x)))
      }
      med <- stats::median(x[k])
      dev <- stats::mad(x[k], center = med)
      .rlna(
        if (dev > 0) abs(x - med) <= r[["n"]] * dev else x == med,
        x,
        napass
      )
    },
    cmpl = if (ncol(v)) stats::complete.cases(v) else TRUE,
    dstn = {
      have <- if (is.null(r[["cols"]])) {
        unique(names(v))
      } else {
        intersect(r[["cols"]], names(v))
      }
      if (length(have)) !duplicated(v[, have, drop = FALSE]) else TRUE
    },
    colcnt = ncol(v) == r[["n"]]
  )
}

#' One rule on a table: NULL when it passes or does not apply, else the
#' issue.
#' @noRd
.rlrun <- function(r, data) {
  col <- r[["col"]]
  if (!is.null(col) && !col %in% names(data)) {
    return(NULL)
  }
  res <- tryCatch(
    .rlans(r, if (is.null(col)) data else data[[col]]),
    error = function(e) e
  )
  iss <- function(m, ...) {
    list(severity = r[["sev"]], message = m, ...)
  }
  if (inherits(res, "error")) {
    return(iss(sprintf(
      "%s (the rule could not be evaluated: %s)",
      r[["msg"]],
      conditionMessage(res)
    )))
  }
  if (is.logical(res) && !length(res) && !nrow(data)) {
    return(NULL)
  }
  if (!is.logical(res) || !length(res)) {
    return(iss(sprintf(
      "%s (the rule returned %s, not TRUE or FALSE)",
      r[["msg"]],
      class(res)[1L]
    )))
  }
  # one answer for the table; on a table of one row, an answer for that
  # row unless the rule is about the table as a whole
  if (
    length(res) == 1L &&
      (nrow(data) != 1L || r[["rule"]] %in% c("incr", "colcnt"))
  ) {
    return(if (isTRUE(res)) NULL else iss(r[["msg"]]))
  }
  bad <- which(is.na(res) | !res)
  if (!length(bad)) {
    return(NULL)
  }
  iss(
    sprintf("%s (%d of %d rows)", r[["msg"]], length(bad), nrow(data)),
    nfail = length(bad),
    rows = utils::head(bad, 10L)
  )
}
