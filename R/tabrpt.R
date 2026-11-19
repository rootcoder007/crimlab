#' One call on a published table: check, change, rates and trend
#'
#' The questions asked of a table an agency publishes are nearly always
#' the same: are its categories sound, what changed from one period to
#' the next and how sure can one be, what are the rates where there is a
#' population, and does a short series trend. `tabrpt()` answers them in
#' one call on a table in long form (one row per period and group), with
#' the steps a reviewer needs and that are easy to forget on a deadline:
#' exact intervals and p-values adjusted over every comparison, and the
#' envelope that rounding in the release allows.
#'
#' The steps: [catchk()] on the group columns (with `by`); [yoypx()] on
#' the counts by period within each group, with `rnd`, `kind`, `level`,
#' `adj` and `alpha`; with `pop`, [grprr()] of the counts over the
#' population pooled by group (by period when there is no `by`), per
#' `per`; and for each group with three periods or more and a count
#' above 0, [trndrr()] of its totals by period (over its population by
#' period with `pop`). A group is the combination of its `by` columns,
#' labelled with them joined by `" / "`. The refusals and warnings of the
#' steps are given under the name `tabrpt`.
#'
#' @param data The table: a data frame in long form.
#' @param value The column of the counts: one name.
#' @param prd The column of the periods (numbers: a year, a month index):
#'   one name.
#' @param by `NULL`, or the columns of the groups: names.
#' @param pop `NULL`, or the column of the populations: one name.
#' @param rnd `NULL`, or the rounding base of the published counts.
#' @inheritParams pbbnd
#' @param per The rates are given per `per` of `pop`: one positive number.
#' @param level The confidence level of every interval.
#' @param adj The adjustment of the p-values for the number of
#'   comparisons, as [stats::p.adjust()] names it.
#' @param alpha The level at which an adjusted p-value is marked.
#' @return An object of class `tabrpt`, a list: `check` (the [catchk()]
#'   result, `NULL` without `by`), `change` (the [yoypx()] table),
#'   `rates` (the [grprr()] result, `NULL` without `pop`) and `trend` (a
#'   list of [trndrr()] results named by group, `NULL` for a group too
#'   short to test).
#' @examples
#' d <- data.frame(year = rep(2019:2023, 2),
#'                 svc = rep(c("A", "B"), each = 5),
#'                 n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33),
#'                 pop = rep(c(50000, 12000), each = 5))
#' r <- tabrpt(d, "n", "year", by = "svc", pop = "pop", rnd = 5)
#' r
#' r$change[, c("group", "prd", "pctchg", "padj")]
#' @export
tabrpt <- function(data, value, prd, by = NULL, pop = NULL, rnd = NULL,
                   kind = "nearest", per = 1000, level = 0.95, adj = "BH",
                   alpha = 0.05) {
  f <- "tabrpt"
  if (!is.data.frame(data) || !nrow(data)) {
    .ckstop(f, "`data` must be a data frame with rows")
  }
  .tbcol(value, "value", data, f)
  .tbcol(prd, "prd", data, f)
  if (!is.null(by)) .tbcol(by, "by", data, f, many = TRUE)
  if (!is.null(pop)) .tbcol(pop, "pop", data, f)
  g <- NULL
  if (!is.null(by)) {
    g <- do.call(paste, c(unname(as.list(data[by])), sep = " / "))
  }
  x <- data[[value]]
  p <- data[[prd]]
  .clas(f, {
    check <- if (is.null(by)) NULL else catchk(data[by])
    change <- yoypx(x, p, by = g, rnd = rnd, kind = kind, level = level,
                    adj = adj, alpha = alpha)
    rates <- if (is.null(pop)) {
      NULL
    } else {
      grprr(x, data[[pop]], grp = if (is.null(g)) as.character(p) else g,
            per = per, level = level)
    }
    trend <- .tbtrnd(x, p, g, if (is.null(pop)) NULL else data[[pop]], level)
  })
  structure(list(check = check, change = change, rates = rates,
                 trend = trend), class = "tabrpt")
}

#' @rdname tabrpt
#' @param x A `tabrpt` result.
#' @inheritParams freqtbl
#' @export
print.tabrpt <- function(x, ...) {
  ch <- x[["change"]]
  cat("tabrpt: ", nrow(ch), " group-periods, ", sum(!is.na(ch[["p"]])),
      " changes tested, ", sum(ch[["signif"]], na.rm = TRUE),
      " marked after adjustment\n", sep = "")
  if (!is.null(x[["check"]])) {
    ok <- isTRUE(x[["check"]][["ok"]])
    cat("groups: ", if (ok) "no problem found" else "see $check", "\n",
        sep = "")
  }
  if (!is.null(x[["rates"]])) {
    r <- x[["rates"]][["rates"]]
    cat("rates per ", format(x[["rates"]][["per"]]), ": ",
        paste0(r[["group"]], " ", format(r[["rate"]], digits = 3),
               collapse = ", "), "\n", sep = "")
  }
  for (k in names(x[["trend"]])) {
    t <- x[["trend"]][[k]]
    if (is.null(t)) {
      cat("trend ", k, ": too short to test\n", sep = "")
    } else {
      cat("trend ", k, ": rate ratio per period ",
          format(t[["rr"]][["rr"]], digits = 3), ", p ",
          format(t[["rr"]][["p"]], digits = 3), "\n", sep = "")
    }
  }
  invisible(x)
}

#' A column argument: one name (many: one or more) of a column of data.
#' @noRd
.tbcol <- function(v, what, data, f, many = FALSE) {
  if (!is.character(v) || anyNA(v) || !length(v) ||
        (!many && length(v) != 1L) || !all(v %in% names(data))) {
    .ckstop(f, "`%s` must name %s of `data`", what,
            if (many) "columns" else "one column")
  }
}

#' The trend of each group (all rows without groups g): NULL when it has
#' fewer than three periods or no count.
#' @noRd
.tbtrnd <- function(x, p, g, pv, level) {
  gg <- if (is.null(g)) rep("all", length(x)) else g
  trend <- lapply(unique(gg), function(k) {
    i <- gg == k
    t <- sort(unique(p[i]))
    y <- vapply(t, function(s) sum(x[i][p[i] == s]), 0)
    e <- if (is.null(pv)) {
      NULL
    } else {
      vapply(t, function(s) sum(pv[i][p[i] == s]), 0)
    }
    if (length(t) < 3L || sum(y) <= 0) {
      return(NULL)
    }
    trndrr(y, t, expo = e, level = level)
  })
  names(trend) <- unique(gg)
  trend
}
