#' Placement volatility of each person-year
#'
#' Counts how often a person's placement region changes within a period
#' (a fiscal year), from records that each carry the region at the time of
#' placement `ra` and the region of the most recent placement `rb`. For a
#' person-year with records `r = 1, ..., R` in a given order the count is
#' \deqn{v = \sum_r 1(A_r \ne B_r) + \sum_{r > 1} 1(A_r \ne A_{r-1}),}{v =
#' sum_r 1(A_r != B_r) + sum_{r > 1} 1(A_r != A_{r-1}),}
#' a change within a record (`within`) plus a change between consecutive
#' records (`across`), as Ruhela (2026) defines placement volatility for
#' Ontario's restrictive confinement records (Appendix H).
#'
#' The first sum does not depend on the order of the records, the second
#' does. A release without placement dates fixes no order, so `ord` says
#' which is used: `"region"` sorts each person-year's records by `ra` (the
#' order of the published analysis code), under which `across` is the
#' number of distinct regions at placement less one and no longer depends
#' on how the rows arrived; `"given"` keeps the rows in the order given.
#' See `proof("plcvol")`: `0 <= v <= 2R - 1`, with equality possible; `v
#' = 0` when every record stays in one region; the count is unchanged by
#' an injective relabelling of the regions and can fall when two regions
#' are merged.
#'
#' Rows missing `id`, `ra` or `rb` are dropped. The work is one ordering
#' of the rows, `O(n log n)` for `n` rows, then one pass.
#'
#' @param ra,rb The region at the time of placement and the region of the
#'   most recent placement of each record: labels (text, factors or
#'   numbers), compared as text.
#' @param id The person-year of each record: one label, or a list of labels
#'   (person and year) combined.
#' @param ord `"region"` or `"given"`, the order of the records within a
#'   person-year.
#' @return A data frame of class `plcvol`, one row per person-year in the
#'   order of first appearance: `id` (the labels, one column each when
#'   `id` is a list), `n` (records), `within`, `across` and `vm` (their
#'   sum).
#' @references Ruhela, V. S. (2026). Alert complexity and placement
#'   volatility in Ontario restrictive confinement data. Major research
#'   paper, Centre for Criminology and Sociolegal Studies, University of
#'   Toronto.
#' @examples
#' ra <- c("East", "East", "West", "North", "North")
#' rb <- c("East", "West", "West", "North", "East")
#' plcvol(ra, rb, id = c(1, 1, 1, 2, 2))
#' plcvol(ra, rb, id = c(1, 1, 1, 2, 2), ord = "given")
#' @export
plcvol <- function(ra, rb, id, ord = c("region", "given")) {
  fn <- "plcvol"
  ord <- .ckmarg("ord", fn)
  g <- .pvgrp(id, length(ra), fn)
  if (length(rb) != length(ra)) {
    .ckstop(fn, "`ra` and `rb` must have one value per record")
  }
  a <- as.character(ra)
  b <- as.character(rb)
  ok <- !is.na(a) & !is.na(b) & !is.na(g$key)
  if (!any(ok)) {
    .ckstop(fn, "no record has `id`, `ra` and `rb` all present")
  }
  k <- factor(g$key[ok], unique(g$key[ok]))
  a <- a[ok]
  b <- b[ok]
  o <- if (ord == "region") order(k, a, method = "radix") else order(k)
  a <- a[o]
  k2 <- k[o]
  same <- c(FALSE, k2[-1L] == k2[-length(k2)])
  acr <- as.integer(same & a != c(NA, a[-length(a)]))
  out <- g$lab[ok, , drop = FALSE][!duplicated(k), , drop = FALSE]
  out$n <- as.vector(table(k))
  out$within <- as.vector(tapply(a != b[o], k2, sum))
  out$across <- as.vector(tapply(acr, k2, sum))
  out$vm <- out$within + out$across
  rownames(out) <- NULL
  structure(out, class = c("plcvol", "data.frame"), ord = ord)
}

#' @export
print.plcvol <- function(x, ...) {
  cat(
    "plcvol: ",
    nrow(x),
    " person-years, records in ",
    attr(x, "ord"),
    " order; mean volatility ",
    format(mean(x$vm), digits = 4),
    " (within ",
    format(mean(x$within), digits = 4),
    ", across ",
    format(mean(x$across), digits = 4),
    ")\n",
    sep = ""
  )
  print(utils::head(as.data.frame(unclass(x)), 10L), row.names = FALSE)
  invisible(x)
}

#' Alert complexity of each person-year
#'
#' Counts the distinct alert configurations among a person-year's records.
#' Each record carries flags (in Ontario's records the mental health,
#' suicide risk and suicide watch alerts); its configuration is the tuple
#' of the flags, `(1, 0, 0)` for a mental health alert alone, coded as the
#' text `"100"`. Alert complexity is the number of distinct configurations
#' the person-year's records span, and `cx`, complexity of two or more, is
#' the treatment of Ruhela (2026, Appendix H).
#'
#' `cx = 1` holds exactly when two records are in different
#' configurations: a record without alerts and one with a single alert
#' give `cx = 1`, while one record with all three alerts gives `cx = 0`.
#' The count measures change of alert state over the records, not alerts
#' held at once (`proof("alrtcx")`). It does not depend on the order of the
#' records or on an injective relabelling of the configurations, never
#' falls when a record is added, and is at most the number of
#' configurations allowed. `keep` restricts the count to the
#' configurations given (the published analysis counts the six seen in
#' the release and ignores any other). The work is `O(n log n)`.
#'
#' @param flags The alert flags of each record: a data frame, list or
#'   matrix with one column per flag; each value 1 / 0, `TRUE` / `FALSE`
#'   or yes / no text (case and surrounding space ignored).
#' @param id The person-year of each record, as for [plcvol()].
#' @param keep `NULL` to count every configuration, or the configurations
#'   to count, as codes such as `"100"` (one digit per flag).
#' @return A data frame of class `alrtcx`, one row per person-year in order
#'   of first appearance: `id`, `n` (records), `ac` (distinct
#'   configurations among those counted), `cx` (1 when `ac >= 2`), and
#'   `anyal` (1 when any record carries an alert).
#' @inherit plcvol references
#' @examples
#' f <- data.frame(mh = c("Yes", "No", "No", "Yes"),
#'                 sr = c("No", "No", "Yes", "Yes"),
#'                 sw = c("No", "No", "No", "Yes"))
#' alrtcx(f, id = c(1, 1, 2, 3))
#' @export
alrtcx <- function(flags, id, keep = NULL) {
  fn <- "alrtcx"
  f <- .acflags(flags, fn)
  g <- .pvgrp(id, nrow(f), fn)
  code <- do.call(paste0, as.data.frame(f))
  code[rowSums(is.na(f)) > 0L] <- NA
  ok <- !is.na(g$key) & !is.na(code)
  if (!any(ok)) {
    .ckstop(fn, "no record has `id` and every flag present")
  }
  .cktxt(keep, "keep", fn)
  if (
    !is.null(keep) &&
      (!is.character(keep) || anyNA(keep) || any(nchar(keep) != ncol(f)))
  ) {
    .ckstop(
      fn,
      "`keep` must be codes of %d digits such as \"%s\"",
      ncol(f),
      strrep("1", ncol(f))
    )
  }
  k <- factor(g$key[ok], unique(g$key[ok]))
  cnt <- if (is.null(keep)) {
    code[ok]
  } else {
    ifelse(code[ok] %in% keep, code[ok], NA)
  }
  out <- g$lab[ok, , drop = FALSE][!duplicated(k), , drop = FALSE]
  out$n <- as.vector(table(k))
  out$ac <- as.vector(tapply(cnt, k, function(v) length(unique(v[!is.na(v)]))))
  out$cx <- as.integer(out$ac >= 2L)
  out$anyal <- as.vector(tapply(rowSums(f[ok, , drop = FALSE]) > 0, k, any)) +
    0L
  rownames(out) <- NULL
  structure(out, class = c("alrtcx", "data.frame"))
}

#' @export
print.alrtcx <- function(x, ...) {
  cat(
    "alrtcx: ",
    nrow(x),
    " person-years; ",
    sum(x$cx),
    " (",
    format(100 * mean(x$cx), digits = 3),
    "%) span two or more ",
    "alert configurations\n",
    sep = ""
  )
  print(utils::head(as.data.frame(unclass(x)), 10L), row.names = FALSE)
  invisible(x)
}

# The person-year of each record: a key (labels joined) and the labels as
# a data frame, from one vector or a list of them, each of length n.
.pvgrp <- function(id, n, fn) {
  l <- if (is.list(id)) as.list(id) else list(id = id)
  if (!length(l) || any(lengths(l) != n) || n == 0L) {
    .ckstop(fn, "`id` must give one label per record (a vector or a list)")
  }
  if (is.null(names(l)) || !all(nzchar(names(l)))) {
    names(l) <- if (length(l) == 1L) "id" else paste0("id", seq_along(l))
  }
  lab <- as.data.frame(l, stringsAsFactors = FALSE, optional = TRUE)
  key <- do.call(.clkey, unname(l))
  key[Reduce(`|`, lapply(l, is.na))] <- NA
  list(key = key, lab = lab)
}

# The flags as a 0/1 integer matrix with NA for a value that is none of
# the accepted forms.
.acflags <- function(flags, fn) {
  if (is.matrix(flags)) {
    flags <- as.data.frame(flags)
  }
  if (
    !is.list(flags) ||
      !length(flags) ||
      length(unique(lengths(flags))) != 1L ||
      !length(flags[[1L]])
  ) {
    .ckstop(fn, "`flags` must be columns of equal length, one per flag")
  }
  .cktxtdf(flags, "flags", fn)
  m <- vapply(flags, .acflag1, integer(length(flags[[1L]])))
  matrix(m, ncol = length(flags))
}

# One column of flags as 0/1 integers, NA for anything else.
.acflag1 <- function(v) {
  if (is.logical(v) || is.numeric(v)) {
    # 0 or 1 exactly: 0.6 is not a flag, and as.integer() would cut it
    out <- rep(NA_integer_, length(v))
    ok <- !is.na(v) & (v == 0 | v == 1)
    out[ok] <- as.integer(v[ok])
    return(out)
  }
  s <- tolower(trimws(as.character(v)))
  # the 1 forms, then the 0 forms; anything else (NA too) is NA
  c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)[
    match(s, c("yes", "y", "true", "1", "no", "n", "false", "0"))
  ]
}
