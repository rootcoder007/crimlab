# Period-over-period change tables out: labels, a summary, and writers to
# CSV, TSV, JSON, Markdown, HTML and PDF. The input is a change table:
# yoypx() gives one, and any data frame with its columns (prd, value,
# prev, change and pctchg or ppchg; group, flag, pctlo and pctup when
# there are any) is one. Numbers are written by the compiled core
# (src/ynum.h): the delimited and JSON writers write every digit, so a
# value read back is the value written; the display writers round only
# what they show.

#' Label the periods of a change table without changing them
#'
#' A period's value and its name are not the same thing: a fiscal year
#' keyed by the year it ends (2023 for the year that runs from 2022
#' into 2023) printed as "2023" names a calendar year the row is not about.
#' `yoylbl()` adds a column `perlbl` that the writers show in place of the
#' period; the comparison still runs on the value, so relabelling cannot
#' change an answer. `fylbl()` names a fiscal year by the years it spans.
#'
#' `labels` is a function of the periods, a vector named by period (a
#' partial mapping: periods it does not name keep their value), or an
#' unnamed vector with one label per distinct period, in their order of
#' appearance.
#'
#' @param x A change table ([yoypx()]).
#' @param labels A function, a named vector or one label per period.
#' @param endyear The year each fiscal year ends, whole numbers (or text
#'   holding them).
#' @param sep The text between the two years.
#' @param short `TRUE` writes the second year with two digits.
#' @return `yoylbl()`: `x` with a column `perlbl`. `fylbl()`: text, `NA`
#'   where `endyear` is missing.
#' @seealso [yoypx()], [yoysmry()], [yoywr()].
#' @examples
#' fylbl(c(2023, 2024))
#' fylbl(2023, sep = "-", short = FALSE)
#' # simulated admissions to custody by fiscal year
#' y <- yoypx(c(120, 131, 98), 2021:2023)
#' yoylbl(y, fylbl)
#' yoylbl(y, c("2023" = "the last year"))
#' @export
yoylbl <- function(x, labels) {
  d <- .ytab(x, "yoylbl")
  p <- d$prd
  if (!is.function(labels) && (!is.atomic(labels) || is.null(labels))) {
    stop("yoylbl: `labels` must be a function of the periods or a vector ",
         "of labels.", call. = FALSE)
  }
  lab <- if (is.function(labels)) {
    tryCatch(as.character(labels(p)), error = function(e) {
      stop("yoylbl: `labels` failed on the periods: ", conditionMessage(e),
           call. = FALSE)
    })
  } else if (!is.null(names(labels))) {
    out <- as.character(p)
    hit <- match(out, names(labels))
    out[!is.na(hit)] <- as.character(labels)[hit[!is.na(hit)]]
    out
  } else {
    u <- unique(p)
    if (length(labels) != length(u)) {
      stop(sprintf(paste0("yoylbl: `labels` has %d entries for %d distinct ",
                          "periods; name them, or pass a function."),
                   length(labels), length(u)), call. = FALSE)
    }
    as.character(labels)[match(p, u)]
  }
  if (length(lab) != nrow(d)) {
    stop("yoylbl: `labels` must give one label per row.", call. = FALSE)
  }
  x$perlbl <- lab
  x
}

#' @rdname yoylbl
#' @export
fylbl <- function(endyear, sep = "/", short = TRUE) {
  sep <- .ytxt1(sep, "sep", "fylbl")
  short <- .yflag(short, "short", "fylbl")
  if (!is.null(endyear) && !is.numeric(endyear) && !is.character(endyear)) {
    stop("fylbl: `endyear` must be numeric or character years.",
         call. = FALSE)
  }
  .cktxt(endyear, "endyear", "fylbl")
  if (is.numeric(endyear) && any(endyear < 1, na.rm = TRUE)) {
    stop("fylbl: `endyear` must be positive years.", call. = FALSE)
  }
  if (is.numeric(endyear) && any(endyear != round(endyear), na.rm = TRUE)) {
    stop("fylbl: `endyear` must be whole years.", call. = FALSE)
  }
  # a year is converted only when it is finite and an integer holds it:
  # what as.integer() gives for NaN or a huge value differs by platform
  num <- if (is.numeric(endyear)) {
    endyear
  } else {
    suppressWarnings(as.numeric(as.character(endyear)))
  }
  fits <- is.finite(num) & abs(num) <= .Machine$integer.max
  y <- rep(NA_integer_, length(num))
  y[fits] <- as.integer(trunc(num[fits]))
  out <- rep(NA_character_, length(y))
  ok <- !is.na(y)
  tail <- if (short) sprintf("%02d", y[ok] %% 100L) else as.character(y[ok])
  out[ok] <- paste0(y[ok] - 1L, sep, tail)
  out
}

#' Summarise a change table
#'
#' One row per series (per group of a grouped table): the first and last
#' period with a value, the values there, the change from the first to
#' the last in percent, the compound rate per step, and how many of the
#' period-to-period changes went up, down or nowhere.
#'
#' The series are put in period order first. `totalpct` and `cagrpct` are
#' `NA` unless the first value is positive (a compound rate out of zero
#' does not exist); `cagrpct` is \eqn{100((last/first)^{1/s} - 1)} with
#' \eqn{s} = (values - 1) / `lag` steps.
#'
#' @inheritParams yoylbl
#' @param lag The periods per step, the `lag` the table was made with.
#' @return A data frame of class `yoysmry`: `group` (for a grouped table),
#'   `from`, `to`, `first`, `last`, `totalpct`, `cagrpct`, `periods`, `up`,
#'   `down` and `flat`.
#' @seealso [yoypx()], [yoywr()].
#' @examples
#' # simulated police-reported counts of two offences
#' d <- data.frame(year = rep(2019:2023, 2),
#'                 offence = rep(c("Assault", "Robbery"), each = 5),
#'                 n = c(812, 845, 790, 901, 954, 120, 131, 98, 105, 140))
#' yoysmry(yoypx(d$n, d$year, d$offence))
#' @export
yoysmry <- function(x, lag = 1) {
  d <- .ytab(x, "yoysmry")
  .ylag(lag, "yoysmry")
  # NROW() of a column: a table subset to no rows can come back from
  # as.data.frame(unclass()) without a row count
  if (!NROW(d[["prd"]])) {
    stop("yoysmry: `x` has no rows: there is nothing to summarise.",
         call. = FALSE)
  }
  chg <- d[[.ychg(d)]]
  g <- if ("group" %in% names(d)) as.character(d$group) else rep("", nrow(d))
  rows <- lapply(unique(g), function(k) {
    i <- which(g == k)
    i <- i[order(d$prd[i], method = "radix")]
    r <- .pcall(23L, as.numeric(d$value[i]), as.numeric(chg[i]), lag)
    fin <- i[!is.na(d$value[i])]
    out <- data.frame(
      from = if (length(fin)) d$prd[fin[1L]] else NA,
      to = if (length(fin)) d$prd[fin[length(fin)]] else NA,
      first = r[1L], last = r[2L], totalpct = r[3L], cagrpct = r[4L],
      periods = as.integer(r[5L]), up = as.integer(r[6L]),
      down = as.integer(r[7L]), flat = as.integer(r[8L])
    )
    if ("group" %in% names(d)) cbind(group = k, out) else out
  })
  structure(do.call(rbind, rows), class = c("yoysmry", "data.frame"))
}

#' @export
print.yoysmry <- function(x, ...) {
  cat("yoysmry:", nrow(x), "series\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

#' Write a change table as delimited text, JSON or Markdown
#'
#' `yoycsv()` and `yoytsv()` write every column of the table, one row per
#' line, numbers with every digit (the shortest text that reads back as
#' the same number), text as it is (a code such as `046` stays `046`),
#' `TRUE` and `FALSE` as such and missing values as `na`; `meta` lines
#' come first as comments. In CSV a field holding a comma, a quote or a
#' line break is quoted, and an empty text field is written as `""` so
#' that it stays apart from a missing value written as empty. A TSV field
#' cannot hold a tab or a line break, so those become spaces. `yoyjs()`
#' writes a JSON object of the `meta` lines, the summary ([yoysmry()]) and
#' the rows, every digit kept. `yoymd()` writes a Markdown table for
#' reading, with the values rounded for display and Markdown's own
#' characters escaped.
#'
#' The delimited and JSON writers change no value: `digits`, when given,
#' rounds the percentage columns (`pctchg`, `ppchg`, `pctlo`, `pctup`) and
#' nothing else, and only when asked.
#'
#' @param x A change table ([yoypx()]); `yoycsv()`, `yoytsv()` and
#'   `yoyjs()` take any data frame.
#' @param file A file path, or `NULL` to return the text.
#' @param digits `NULL`, or the decimals of the percentage columns (0 to
#'   15); in `yoymd()` the decimals of the change shown.
#' @param na The text of a missing value.
#' @param meta `NULL`, or lines that describe the table (its source, its
#'   settings), written as comments or as the table's caption.
#' @param comment The comment mark of the `meta` line.
#' @param pretty `TRUE` indents the JSON.
#' @param align `TRUE` pads the Markdown columns to a common width.
#' @return The text when `file` is `NULL`, otherwise `file`, invisibly.
#' @seealso [yoywr()], [yoyhtml()].
#' @examples
#' # simulated police-reported counts by offence and year, any agency
#' d <- data.frame(year = rep(2021:2023, 2),
#'                 offence = rep(c("Assault", "Robbery"), each = 3),
#'                 n = c(812, 845, 790, 120, 131, 98))
#' y <- yoypx(d$n, d$year, d$offence)
#' cat(yoycsv(y, meta = "simulated counts"))
#' cat(yoymd(y))
#' # simulated admissions to custody, with a code that must stay "046"
#' a <- data.frame(code = c("046", "112"), n = c(1 / 3, 2))
#' cat(yoytsv(a))
#' cat(yoyjs(y, pretty = FALSE))
#' @name yoycsv
#' @export
yoycsv <- function(x, file = NULL, digits = NULL, na = "", meta = NULL,
                   comment = "#") {
  .ydelim(x, file, ",", digits, na, meta, comment, "yoycsv")
}

#' @rdname yoycsv
#' @export
yoytsv <- function(x, file = NULL, digits = NULL, na = "", meta = NULL,
                   comment = "#") {
  .ydelim(x, file, "\t", digits, na, meta, comment, "yoytsv")
}

#' @rdname yoycsv
#' @export
yoyjs <- function(x, file = NULL, digits = NULL, pretty = TRUE, meta = NULL) {
  file <- .yfile(file, "yoyjs")
  pretty <- .yflag(pretty, "pretty", "yoyjs")
  d <- .yexp(x, digits, "yoyjs")
  m <- .ymeta(meta, "yoyjs")
  ok <- all(c("prd", "value", "prev", "change") %in% names(d)) &&
    any(c("pctchg", "ppchg") %in% names(d))
  # a table of no rows has no summary (yoysmry() refuses it)
  sm <- if (ok && length(d[["prd"]])) as.data.frame(unclass(yoysmry(d)))
  ind <- if (pretty) 2L else NULL
  parts <- c(
    if (length(m)) c("meta", .yjs1(m)),
    if (ok) c("summary", .yjsrows(sm, ind)),
    "rows", .yjsrows(d, ind)
  )
  k <- parts[c(TRUE, FALSE)]
  v <- parts[c(FALSE, TRUE)]
  txt <- if (pretty) {
    paste0("{\n", paste0("  ", vapply(k, .yjs1, ""), ": ", v,
                         collapse = ",\n"), "\n}")
  } else {
    paste0("{", paste0(vapply(k, .yjs1, ""), ":", v, collapse = ","), "}")
  }
  .yout(paste0(txt, "\n"), file)
}

#' @rdname yoycsv
#' @export
yoymd <- function(x, file = NULL, digits = 1L, align = TRUE, meta = NULL) {
  file <- .yfile(file, "yoymd")
  d <- .ytab(x, "yoymd")
  digits <- .ydig(digits, "yoymd")
  align <- .yflag(align, "align", "yoymd")
  m <- .ymeta(meta, "yoymd")
  s <- .yshow(d, digits)
  cols <- s$cols
  cols$note <- ifelse(is.na(s$flag), "", s$flag)
  cols <- lapply(cols, .ymdesc)
  nms <- .ymdesc(names(cols))
  w <- if (align) {
    vapply(seq_along(cols), function(i) {
      max(3L, nchar(c(nms[i], cols[[i]]), type = "width"))
    }, 0L)
  } else {
    rep(0L, length(cols))
  }
  nby <- as.integer("group" %in% names(d))
  right <- seq_along(cols) > nby + 1L & names(cols) != "note"
  pad <- function(s, k) {
    gap <- strrep(" ", pmax(0L, w[k] - nchar(s, type = "width")))
    if (right[k]) paste0(gap, s) else paste0(s, gap)
  }
  row <- function(vals) {
    paste0("| ", paste(vapply(seq_along(vals), function(k) pad(vals[k], k), ""),
                       collapse = " | "), " |")
  }
  rule <- paste0("| ", paste(vapply(seq_along(w), function(k) {
    n <- max(3L, w[k])
    if (right[k]) paste0(strrep("-", n - 1L), ":") else strrep("-", n)
  }, ""), collapse = " | "), " |")
  lines <- c(row(nms), rule, vapply(seq_len(nrow(d)), function(r) {
    row(vapply(cols, function(v) v[r], ""))
  }, ""))
  if (length(m)) {
    lines <- c(lines, "", paste0("_", .ymdesc(m), "_"))
  }
  .yout(paste0(paste(lines, collapse = "\n"), "\n"), file)
}

.yfmt <- function(x, mode, d = 0L) {
  .Call(C_clab_r_fmt, as.numeric(x), as.integer(mode), as.integer(d))
}

.ytab <- function(x, fn) {
  need <- c("prd", "value", "prev", "change")
  if (!is.data.frame(x) || !all(need %in% names(x)) ||
        !any(c("pctchg", "ppchg") %in% names(x))) {
    stop(fn, ": `x` must be a change table: a data frame with the columns ",
         "prd, value, prev, change and pctchg or ppchg, as yoypx() gives.",
         call. = FALSE)
  }
  .cktxtdf(x, "x", fn)
  if (anyDuplicated(names(x))) {
    .ckstop(fn, "`x` has two columns named %s",
            names(x)[anyDuplicated(names(x))])
  }
  as.data.frame(unclass(x))
}

.ychg <- function(d) if ("ppchg" %in% names(d)) "ppchg" else "pctchg"

.ydig <- function(d, fn, nm = "digits", null = FALSE) {
  if (null && is.null(d)) {
    return(NULL)
  }
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || d < 0 || d > 15 ||
        d != floor(d)) {
    stop(fn, ": `", nm, "` must be one whole number from 0 to 15.",
         call. = FALSE)
  }
  as.integer(d)
}

.yflag <- function(x, nm, fn) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(fn, ": `", nm, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

# The lag of a summary: one positive number.
.ylag <- function(lag, fn) {
  if (!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag <= 0) {
    stop(fn, ": `lag` must be one positive number.", call. = FALSE)
  }
  lag
}

.ytxt1 <- function(x, nm, fn) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(fn, ": `", nm, "` must be one string.", call. = FALSE)
  }
  .cktxt(x, nm, fn)
}

.yfile <- function(file, fn) {
  if (!is.null(file) && (!is.character(file) || length(file) != 1L ||
                           is.na(file) || !nzchar(file))) {
    stop(fn, ": `file` must be NULL or one file path.", call. = FALSE)
  }
  .ckpath(file, fn)
}

.yout <- function(txt, file) {
  if (is.null(file)) {
    return(txt)
  }
  writeBin(charToRaw(enc2utf8(txt)), file)
  invisible(file)
}

.ymeta <- function(meta, fn) {
  if (!is.null(meta) && (!is.character(meta) || anyNA(meta))) {
    stop(fn, ": `meta` must be NULL or text.", call. = FALSE)
  }
  .cktxt(meta, "meta", fn)
  if (length(meta)) gsub("[\r\n]+", " ", paste(meta, collapse = "; "))
}

# the period as shown: its label (yoylbl()) when there is one
.yper <- function(d) {
  as.character(if ("perlbl" %in% names(d)) d$perlbl else d$prd)
}

# decimals of the values shown: none for counts (a table with an exact
# interval), else the fewest that show every value to six places
.yvd <- function(d) {
  if ("pctlo" %in% names(d)) {
    return(0L)
  }
  v <- d$value[is.finite(d$value)]
  if (!length(v)) {
    return(0L)
  }
  min(6L, max(vapply(v, function(z) {
    for (k in 0:6) if (abs(z - round(z, k)) < 1e-9) return(k)
    6L
  }, 0L)))
}

.ynum <- function(v, dg) {
  out <- .yfmt(v, 1L, dg)
  out[is.na(v)] <- "\u2014"
  out
}

.ypct <- function(v, dg, pp) {
  out <- .yfmt(v, if (pp) 3L else 2L, dg)
  out[is.na(v) | is.nan(v)] <- "\u2014"
  out
}

# control characters shown escaped, never written; ws keeps tabs and line
# breaks
.yctl <- function(s, ws = FALSE) {
  s <- as.character(s)
  u <- function(k) intToUtf8(k)
  re <- paste0("[\001-\037\177", u(0x80), "-", u(0x9f), u(0x200e),
               u(0x200f), u(0x202a), "-", u(0x202e), u(0x2066), "-",
               u(0x2069), "]")
  hit <- which(!is.na(s) & grepl(re, s, perl = TRUE))
  nm <- c(`7` = "\\a", `8` = "\\b", `9` = "\\t", `10` = "\\n", `12` = "\\f",
          `13` = "\\r", `27` = "\\033")
  for (i in hit) {
    cp <- utf8ToInt(enc2utf8(s[i]))
    if (anyNA(cp)) next
    s[i] <- paste(vapply(cp, function(k) {
      if (ws && k %in% c(9L, 10L)) {
        intToUtf8(k)
      } else if (as.character(k) %in% names(nm)) {
        nm[[as.character(k)]]
      } else if (k < 32L || k == 127L) {
        sprintf("\\%03o", k)
      } else if (k >= 0x80 && k <= 0x9f || k %in% c(0x200e, 0x200f) ||
                   k >= 0x202a && k <= 0x202e || k >= 0x2066 && k <= 0x2069) {
        sprintf("<U+%04X>", k)
      } else {
        intToUtf8(k)
      }
    }, ""), collapse = "")
  }
  s
}

# a data frame with the percentage columns rounded when digits is given
.yexp <- function(x, digits, fn) {
  if (!is.data.frame(x)) {
    stop(fn, ": `x` must be a data frame.", call. = FALSE)
  }
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  digits <- .ydig(digits, fn, null = TRUE)
  if (!is.null(digits)) {
    for (nm in intersect(c("pctchg", "ppchg", "pctlo", "pctup"), names(d))) {
      d[[nm]] <- round(d[[nm]], digits)
    }
  }
  d
}

# one column as text: numbers with every digit
.ycol <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.numeric(v)) .yfmt(v, 0L) else as.character(v)
}

.ydelim <- function(x, file, sep, digits, na, meta, comment, fn) {
  file <- .yfile(file, fn)
  na <- .ytxt1(na, "na", fn)
  comment <- .ytxt1(comment, "comment", fn)
  d <- .yexp(x, digits, fn)
  m <- .ymeta(meta, fn)
  cell <- function(v, txt) {
    out <- ifelse(is.na(v), na, v)
    if (sep != ",") {
      return(gsub("[\t\r\n]", " ", out))
    }
    q <- !is.na(v) & (grepl("[,\"\r\n]", v) | (txt & !nzchar(v)))
    out[q] <- paste0("\"", gsub("\"", "\"\"", v[q], fixed = TRUE), "\"")
    out
  }
  cols <- lapply(d, function(v) cell(.ycol(v), is.character(v) || is.factor(v)))
  lines <- c(
    if (length(m)) paste0(comment, " ", m),
    paste(cell(names(d), TRUE), collapse = sep),
    if (nrow(d)) do.call(paste, c(unname(cols), sep = sep))
  )
  .yout(paste0(paste(lines, collapse = "\n"), "\n"), file)
}

# JSON text of one value
.yjs1 <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.null(v) || length(v) == 0L || is.na(v)) {
    return("null")
  }
  if (is.logical(v)) {
    return(if (v) "true" else "false")
  }
  if (is.numeric(v)) {
    if (!is.finite(v)) {
      return(paste0("\"", .yfmt(v, 0L), "\""))
    }
    return(.yfmt(v, 0L))
  }
  s <- enc2utf8(as.character(v))
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  cp <- utf8ToInt(s)
  if (any(cp < 32L)) {
    s <- paste(vapply(cp, function(k) {
      if (k < 32L) sprintf("\\u%04x", k) else intToUtf8(k)
    }, ""), collapse = "")
  }
  paste0("\"", s, "\"")
}

# an array of row objects; ind the indentation of the array, NULL flat
.yjsrows <- function(d, ind) {
  if (!length(d) || !NROW(d[[1L]])) {
    return("[]")
  }
  rows <- vapply(seq_len(nrow(d)), function(i) {
    kv <- vapply(names(d), function(k) {
      paste0(.yjs1(k), if (is.null(ind)) ":" else ": ", .yjs1(d[[k]][i]))
    }, "")
    if (is.null(ind)) {
      paste0("{", paste(kv, collapse = ","), "}")
    } else {
      p <- strrep(" ", ind + 4L)
      paste0(strrep(" ", ind + 2L), "{\n", p,
             paste(kv, collapse = paste0(",\n", p)), "\n",
             strrep(" ", ind + 2L), "}")
    }
  }, "")
  if (is.null(ind)) {
    paste0("[", paste(rows, collapse = ","), "]")
  } else {
    paste0("[\n", paste(rows, collapse = ",\n"), "\n", strrep(" ", ind), "]")
  }
}

.ymdesc <- function(v) {
  v <- .yctl(gsub("([\\\\`*])", "\\\\\\1", gsub("[\r\n]+", " ", v)))
  v <- gsub("&", "&amp;", v, fixed = TRUE)
  v <- gsub("<", "&lt;", v, fixed = TRUE)
  v <- gsub(">", "&gt;", v, fixed = TRUE)
  gsub("|", "\\|", v, fixed = TRUE)
}

# the columns shown by the display writers, as text
.yshow <- function(d, digits, intv = TRUE, md = TRUE) {
  vd <- .yvd(d)
  chg <- .ychg(d)
  pp <- chg == "ppchg"
  cols <- list()
  if ("group" %in% names(d)) cols$group <- as.character(d$group)
  cols$prd <- .yper(d)
  cols$value <- .ynum(d$value, vd)
  cols$previous <- .ynum(d$prev, vd)
  cols$change <- .ynum(d$change, vd)
  cols[[if (pp) "points" else "change %"]] <- .ypct(d[[chg]], digits, pp)
  if (intv && all(c("pctlo", "pctup") %in% names(d))) {
    cols$interval <- ifelse(is.na(d$pctlo), if (md) "" else "\u2014",
                            sprintf("[%s, %s]", .ypct(d$pctlo, 0L, FALSE),
                                    .ypct(d$pctup, 0L, FALSE)))
  }
  fl <- if ("flag" %in% names(d)) d$flag else rep(NA_character_, nrow(d))
  list(cols = cols, flag = fl, pp = pp)
}
