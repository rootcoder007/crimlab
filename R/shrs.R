#' Share of a total, in percent
#'
#' Expresses each count as a percentage of a total (for example each
#' division's share of all stops, or each offence's share of all
#' admissions), with a Wilson confidence interval. A share's denominator is
#' the total of the same events, so shares over a complete grouping sum
#' to 100. Use a rate (such as [grprr()]) instead when the denominator is a
#' population; "40% of stops" and "40 stops per 1,000 residents" are
#' different claims.
#'
#' For count \eqn{y} out of total \eqn{n}, the share is \eqn{100\,y/n} and
#' the interval is Wilson's score interval for \eqn{\hat p = y/n},
#' \deqn{\frac{\hat p + \frac{z^2}{2n} \pm z\sqrt{\frac{\hat p(1-\hat
#'   p)}{n} + \frac{z^2}{4n^2}}}{1 + z^2/n},}
#' with \eqn{z} the standard normal quantile for `conflvl`, in percent and
#' clipped to \eqn{[0, 100]} (exactly 0 at a count of 0 and 100 at the
#' whole total). The normal interval on a proportion is wrong in exactly
#' the cases people reach for it (small counts, shares near 0 or 100),
#' where it runs past the ends of the scale; Wilson's stays inside and is
#' accurate there (Wilson, 1927). It treats each count as a binomial draw
#' from the total, an approximation when the total is the sum of the
#' counts.
#'
#' With a data frame, counts are summed within each combination of the
#' `by` columns (missing values ignored), the groups sorted as
#' [stats::aggregate()] sorts them. The default method treats each element
#' of `x` as a row; a named `x` gives a leading `label` column.
#'
#' Errors: a count column that is not numeric or not found; negative
#' counts; a total of zero or below; a count larger than the total;
#' `conflvl` outside (0, 1).
#'
#' @param x A data frame, or a vector of counts.
#' @param count The column of counts (required for a data frame), unquoted
#'   or as a string: non-negative numbers.
#' @param by `NULL`, or the grouping columns (names).
#' @param total The denominator: `NULL` (the sum of the counts, so the
#'   shares sum to 100), one non-negative number, or (for a data frame) the
#'   name of a column whose sum is the total, for when some of it is not in
#'   the table.
#' @param conflvl The level of the interval, in (0, 1).
#' @param ... Unused; anything given is an error.
#' @return A data frame of class `shrs`: the `by` columns (or `label`),
#'   `count`, `total`, `share`, `lower` and `upper` (percent), with an
#'   attribute `share`, a list of `count`, `by`, `total` and `conflvl`.
#' @references Wilson EB (1927). Probable inference, the law of succession,
#'   and statistical inference. *Journal of the American Statistical
#'   Association* 22(158), 209--212. \doi{10.1080/01621459.1927.10502953}
#' @examples
#' # simulated stops by division
#' stops <- data.frame(division = c("North", "South", "East"),
#'                     stops = c(412, 77, 3))
#' shrs(stops, stops, by = "division")
#' shrs(stops, "stops", by = "division", total = 10000)
#' shrs(c(murder = 3, robbery = 7))
#' try(shrs(c(5, 2), total = 4))
#' @export
shrs <- function(x, ...) UseMethod("shrs")

#' @rdname shrs
#' @export
shrs.data.frame <- function(x, count = NULL, by = NULL, total = NULL,
                            conflvl = 0.95, ...) {
  .shnodots(...)
  sub <- substitute(count)
  nm <- if (is.name(sub) && as.character(sub) %in% names(x)) {
    as.character(sub)
  } else {
    tryCatch(count, error = function(e) NULL)
  }
  .shnm(nm, x)
  if (!is.null(by)) by <- .shby(by, x)
  agg <- .shagg(x, nm, by)
  den <- .shden(total, x, agg)
  .shout(agg, den, conflvl, nm, by)
}

#' @rdname shrs
#' @export
shrs.default <- function(x, total = NULL, conflvl = 0.95, ...) {
  .shnodots(...)
  if (!is.atomic(x) || !is.numeric(x) || all(is.na(x)) || any(is.infinite(x))) {
    stop("shrs: `x` must be finite counts (a numeric vector) or a data ",
         "frame of them.", call. = FALSE)
  }
  lab <- names(x)
  agg <- data.frame(count = as.numeric(x))
  if (!is.null(lab)) {
    agg <- cbind(data.frame(label = lab), agg)
  }
  den <- if (is.null(total)) sum(x, na.rm = TRUE) else .shtot(total, "")
  .shout(agg, den, conflvl, "count", NULL)
}

#' @export
print.shrs <- function(x, ...) {
  m <- attr(x, "share")
  cat("shrs: shares of ", format(m$total), " in percent, ",
      format(100 * m$conflvl), "% Wilson intervals\n", sep = "")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

# shrs.data.frame(): `count` names one numeric column of `x`
.shnm <- function(nm, x) {
  if (!is.character(nm) || length(nm) != 1L || is.na(nm)) {
    stop("shrs: `count` must name the column of counts.", call. = FALSE)
  }
  if (!nm %in% names(x)) {
    stop("shrs: `count` column not found: ", nm, ".", call. = FALSE)
  }
  if (!is.numeric(x[[nm]])) {
    stop("shrs: `", nm, "` must be numeric.", call. = FALSE)
  }
}

.shby <- function(by, x) {
  by <- as.character(by)
  miss <- setdiff(by, names(x))
  if (length(miss)) {
    stop("shrs: `by` column not found: ", toString(miss), ".", call. = FALSE)
  }
  by
}

# shrs.data.frame(): the counts summed, by the `by` columns
.shagg <- function(x, nm, by) {
  if (is.null(by)) {
    data.frame(count = sum(x[[nm]], na.rm = TRUE))
  } else {
    a <- stats::aggregate(x[nm], by = x[by], FUN = sum, na.rm = TRUE)
    names(a) <- c(by, "count")
    a
  }
}

# shrs.data.frame(): the denominator: the counts' sum, a column's sum or a
# number
.shden <- function(total, x, agg) {
  if (is.null(total)) {
    sum(agg$count)
  } else if (is.character(total) && length(total) == 1L && !is.na(total)) {
    if (!total %in% names(x)) {
      stop("shrs: `total` column not found: ", total, ".", call. = FALSE)
    }
    sum(x[[total]], na.rm = TRUE)
  } else {
    .shtot(total, "a column name or ")
  }
}

.shnodots <- function(...) {
  if (...length()) {
    stop("shrs: unused argument(s): ", toString(names(list(...))), ".",
         call. = FALSE)
  }
}

.shtot <- function(total, also) {
  if (!is.numeric(total) || length(total) != 1L || !is.finite(total) ||
        total < 0) {
    stop("shrs: `total` must be NULL, ", also,
         "one finite non-negative number.", call. = FALSE)
  }
  as.numeric(total)
}

.shout <- function(agg, den, conflvl, nm, by) {
  if (!is.numeric(conflvl) || length(conflvl) != 1L || is.na(conflvl) ||
        conflvl <= 0 || conflvl >= 1) {
    stop("shrs: `conflvl` must be one number strictly between 0 and 1.",
         call. = FALSE)
  }
  cnt <- agg$count
  if (any(cnt < 0, na.rm = TRUE)) {
    stop("shrs: `count` must be non-negative.", call. = FALSE)
  }
  if (is.na(den) || den <= 0) {
    stop("shrs: the total must be positive.", call. = FALSE)
  }
  if (any(cnt > den, na.rm = TRUE)) {
    stop("shrs: a count is larger than the total, so these are not shares ",
         "of it.", call. = FALSE)
  }
  ok <- !is.na(cnt)
  r <- matrix(NA_real_, length(cnt), 3L)
  r[ok, ] <- .pcall(18L, cnt[ok], a = c(den, conflvl))
  out <- agg[setdiff(names(agg), "count")]
  out$count <- cnt
  out$total <- den
  out$share <- r[, 1L]
  out$lower <- r[, 2L]
  out$upper <- r[, 3L]
  rownames(out) <- NULL
  structure(out, class = c("shrs", "data.frame"),
            share = list(count = nm, by = by, total = den, conflvl = conflvl))
}
