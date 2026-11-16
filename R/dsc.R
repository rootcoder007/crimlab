# Data description: the look at a table before trusting it -- the shape of
# every column, where the gaps are, which values dominate, which rows do
# not belong. Counting and arithmetic are in the compiled core (dsc.cpp);
# data frames are taken apart and put together here.

#' Inline histogram for a numeric vector
#'
#' Draws a one-line sketch of a column's distribution, as block characters
#' (or ASCII where the console cannot render them). The point is density of
#' information: a mean and a standard deviation cannot tell you a column is
#' bimodal, skewed or has a gap, and this can, in the width of a table cell.
#' \code{\link{profcols}} shows one for every numeric column.
#'
#' @details
#' Non-finite values (\code{NA}, \code{NaN}, \code{Inf}) are dropped. The
#' range of the rest is cut into \code{bins} equal-width bins (as
#' \code{graphics::hist()} with those breaks) and each bin's count is drawn
#' as one of eight Unicode block heights, scaled to the largest bin; an empty
#' bin is a space, and a bin with any values is at least the lowest block.
#' When the session is not UTF-8, or the environment variable
#' \code{CRIMLAB_ASCII_ONLY} is set, five ASCII levels (\code{. : - = #}) are
#' used instead. A column with no spread gives one full block in the middle;
#' no finite values give spaces.
#'
#' Errors: \code{x} that is not numeric; \code{bins} that is not one positive
#' whole number.
#'
#' @param x Numeric vector drawn as an inline histogram.
#' @param bins Number of bins, and so of characters, one positive whole
#'   number (default 10), at most 67108864 (2^26).
#' @return A length-1 character vector of \code{bins} characters.
#' @seealso \code{\link{profcols}}, which includes one per numeric column.
#' @examples
#' set.seed(1)
#' # A symmetric distribution peaks in the middle (simulated values).
#' inlhist(stats::rnorm(1000))
#'
#' # A skewed one leans left.
#' inlhist(stats::rexp(1000))
#'
#' # Bimodality is visible here and in no single summary number.
#' inlhist(c(stats::rnorm(500, -3), stats::rnorm(500, 3)))
#'
#' # A constant column has no spread to show.
#' inlhist(rep(5, 10))
#'
#' # more bins for more detail
#' inlhist(stats::rexp(1000), bins = 30)
#'
#' # text is not a distribution
#' try(inlhist(c("a", "b")))
#' @export
inlhist <- function(x, bins = 10L) {
  f <- "inlhist"
  if (is.null(x)) {
    x <- numeric(0)
  } else if (is.factor(x) || is.character(x) || is.list(x) ||
               !(is.numeric(x) || is.logical(x))) {
    .ckstop(f, "`x` must be numeric, got %s", class(x)[1L])
  }
  bins <- .ckpos(bins, "bins", f)
  if (bins > 2^26) .ckstop(f, "`bins` must be at most 67108864")
  .spark(as.numeric(x), bins)
}

#' Frequency table for one column
#'
#' Counts how often each value of a column occurs, with percentages of all
#' rows and of the non-missing rows, and missing values counted as their own
#' row. It is the first look at a categorical variable (a charge type, a
#' region, a response category).
#'
#' @details
#' Two percentage columns, because both questions get asked and conflating
#' them is how missingness gets hidden: \code{pct} is the share of ALL rows,
#' \code{pctvalid} the share of rows where the value is present. When a column
#' is 40\% missing those two differ enormously, and only the second describes
#' the values that are actually there.
#'
#' Values are tabulated as text (\code{as.character()}), so numbers and
#' factors are counted by their printed value, and factor levels that do not
#' occur get no row. Rows are sorted by descending count (ties in the order
#' \code{table()} gives) or kept in sorted value order with
#' \code{sort = FALSE}. Beyond \code{maxlvl} distinct values, the rest are
#' folded into one \code{"(other)"} row so the percentages still sum to 100.
#' The missing row, when there is one, comes last, with \code{value}
#' \code{NA} and \code{pctvalid} \code{NA}. An empty input gives a table with
#' no rows. The print method adds a bar column.
#'
#' Errors: a data frame without \code{col}; a \code{col} that is not in
#' \code{data}; \code{maxlvl} not one positive whole number.
#'
#' @param data A data frame, or a vector.
#' @param col Column name (a string), when \code{data} is a data frame;
#'   ignored for a vector.
#' @param maxlvl Maximum value rows to return, most frequent first (default
#'   25). The remainder are folded into one \code{"(other)"} row so the
#'   percentages still sum to 100.
#' @param sort Logical: sort by descending count (default \code{TRUE});
#'   \code{FALSE} keeps the sorted order of the values.
#' @param x A result of \code{freqtbl()}.
#' @param ... Not used.
#' @return A data frame of class \code{c("freqtbl", "data.frame")} with
#'   columns \code{value} (character), \code{n} (integer count), \code{pct}
#'   (numeric, percent of all rows) and \code{pctvalid} (numeric, percent of
#'   non-missing rows; \code{NA} on the missing row).
#' @seealso \code{\link{profcols}} for every column at once.
#' @examples
#' df <- data.frame(grade = c("a", "b", "b", "c", NA, "b"))
#' freqtbl(df, "grade")
#'
#' # The two percentage columns differ exactly by the missingness.
#' freqtbl(df, "grade")[, c("pct", "pctvalid")]
#'
#' # A vector works directly.
#' freqtbl(c(1, 1, 2, 3, 3, 3))
#'
#' # Sorted value order rather than frequency order.
#' freqtbl(c("c", "a", "b", "a"), sort = FALSE)
#'
#' # Long tails are folded so the percentages still total 100.
#' set.seed(1)
#' ft <- freqtbl(sample(letters, 500, TRUE), maxlvl = 5)
#' ft
#' sum(ft$pct)
#'
#' # a data frame needs a column name
#' try(freqtbl(df))
#' @export
freqtbl <- function(data, col = NULL, maxlvl = 25L, sort = TRUE) {
  f <- "freqtbl"
  .cknofn(data, "data", f)
  .cknofn(col, "col", f)
  .ckflag(sort, "sort", f)
  v <- if (is.data.frame(data)) {
    if (is.null(col)) {
      .ckstop(f, "`col` is required when `data` is a data frame")
    }
    .colnm(col, "col", f)
    if (!col %in% names(data)) .ckstop(f, "no such column: %s", col)
    data[[col]]
  } else {
    .cktxt(data, "data", f)
  }
  maxlvl <- .ckpos(maxlvl, "maxlvl", f)
  out <- data.frame(value = character(0), n = integer(0), pct = numeric(0),
                    pctvalid = numeric(0), stringsAsFactors = FALSE)
  if (length(v)) {
    s <- as.character(v[!is.na(v)])
    lv <- sort(unique(s))
    r <- .Call(C_clab_r_freq, match(s, lv), length(lv),
               as.integer(sum(is.na(v))), as.numeric(length(v)),
               isTRUE(sort), maxlvl, TRUE)
    vals <- c(lv, "(other)")[r[[1L]]]
    vals[is.na(r[[1L]])] <- NA_character_
    out <- data.frame(value = vals, n = r[[2L]], pct = r[[3L]],
                      pctvalid = r[[4L]], stringsAsFactors = FALSE)
  }
  class(out) <- c("freqtbl", "data.frame")
  out
}

#' @rdname freqtbl
#' @export
print.freqtbl <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  # a subset of the columns (x[, c("pct", "pctvalid")]) prints as it is
  if (is.null(d$n)) {
    print(d, row.names = FALSE)
    return(invisible(x))
  }
  if (nrow(d)) {
    d$bar <- strrep("#", as.integer(round(20 * d$n / max(c(1, d$n)))))
  }
  cat("Frequency table:", sum(x$n), "values\n")
  print(d, row.names = FALSE)
  invisible(x)
}

#' Full pairwise correlation table
#'
#' Computes the correlation of every pair of numeric columns and returns it
#' in long form -- one row per pair, which is easier to sort, filter and join
#' than a matrix -- with the number of complete pairs each figure rests on.
#'
#' @details
#' Non-numeric columns are ignored. For each pair, rows where either value is
#' missing are dropped (pairwise deletion), so different pairs can rest on
#' different rows; \code{nprs} says how many. Spearman's rank correlation (the
#' default) is the Pearson correlation of the ranks, robust to outliers and
#' to monotone transformations; see \code{\link{topcors}} for why it is the
#' default. A pair with fewer than \code{minprs} complete rows, or with a
#' constant column (no variance), gets \code{NA}. The pairs are listed in
#' column order (\code{utils::combn()}).
#'
#' Errors: \code{data} not a data frame; fewer than two numeric columns;
#' \code{minprs} not a whole number of at least 2; an invalid \code{meth}.
#'
#' @param data A data frame; its numeric columns are correlated in pairs.
#' @param meth \code{"spearman"} (default) or \code{"pearson"}.
#' @param minprs Minimum complete pairs required before a correlation is
#'   computed, a whole number of at least 2 (default 3). Below it the pair is
#'   \code{NA} rather than a number computed from almost nothing.
#' @param x A result of \code{cortbl()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("cortbl", "data.frame")} with one
#'   row per pair of numeric columns and the columns \code{x} and \code{y}
#'   (character, the column names), \code{correlation} (numeric, in
#'   \eqn{[-1, 1]} or \code{NA}) and \code{nprs} (integer); the method is in
#'   \code{attr(, "method")}.
#' @seealso \code{\link{topcors}} for just the strongest.
#' @examples
#' set.seed(1)
#' df <- data.frame(a = stats::rnorm(100), b = stats::rnorm(100))
#' df$c <- df$a + stats::rnorm(100, sd = 0.2) # simulated
#'
#' cortbl(df)
#'
#' # Pearson instead of Spearman
#' cortbl(df, meth = "pearson")
#'
#' # nprs shows how much data each figure rests on.
#' gappy <- df
#' gappy$a[1:80] <- NA
#' cortbl(gappy)
#'
#' # A pair with too little overlap is NA, not a number from nothing.
#' thin <- df
#' thin$a[1:99] <- NA
#' cortbl(thin)$correlation
#'
#' # one numeric column cannot be correlated
#' try(cortbl(data.frame(a = 1:5, b = letters[1:5])))
#' @export
cortbl <- function(data, meth = c("spearman", "pearson"), minprs = 3L) {
  f <- "cortbl"
  meth <- .ckarg(meth, c("spearman", "pearson"), "meth", f)
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  num <- names(data)[vapply(data, is.numeric, logical(1))]
  if (length(num) < 2L) {
    .ckstop(f, "need at least two numeric columns to correlate")
  }
  minprs <- .ckn(minprs, "minprs", f, integer = TRUE)
  if (length(minprs) != 1L || is.na(minprs) || minprs < 2L) {
    .ckstop(f, "`minprs` must be a single integer of at least 2")
  }
  r <- .Call(C_clab_r_pcor, .nummat(data, num), meth == "pearson", minprs)
  pr <- utils::combn(num, 2L)
  out <- data.frame(x = pr[1L, ], y = pr[2L, ], correlation = r[[1L]],
                    nprs = r[[2L]], stringsAsFactors = FALSE)
  attr(out, "method") <- meth
  class(out) <- c("cortbl", "data.frame")
  out
}

#' @rdname cortbl
#' @export
print.cortbl <- function(x, ...) {
  cat("Correlations (", attr(x, "method"), "): ", nrow(x), " pairs\n",
      sep = "")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(d, "method") <- NULL
  print(d, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Strongest pairwise correlations in a data frame
#'
#' Ranks the numeric column pairs by the strength of their association and
#' returns the strongest, so a wide table's structure can be read without
#' squinting at a correlation matrix. Use it on a new data set to find
#' near-duplicate columns, derived columns and the relationships worth a
#' closer look.
#'
#' @details
#' Spearman is the default deliberately. Pearson measures LINEAR association
#' only, so it understates a relationship that is perfectly monotone but
#' curved, and a single outlier can manufacture or destroy it. On data you
#' have not yet inspected -- which is the situation this function is for --
#' the rank correlation is the safer question to ask.
#'
#' Non-numeric columns are ignored. Each pair uses the rows where both values
#' are present; a pair with fewer than three such rows, or with a constant
#' column, has no correlation and is left out. Pairs are ranked by absolute
#' correlation, those below \code{minabs} dropped, and the first \code{n}
#' returned. \code{\link{cortbl}} gives every pair, with the number of rows
#' behind each.
#'
#' Errors: \code{data} not a data frame; fewer than two numeric columns;
#' \code{n} not one positive whole number; \code{minabs} not one number in
#' \eqn{[0, 1]}; an invalid \code{meth}.
#'
#' @param data A data frame; its numeric pairs are ranked by correlation.
#' @param n Number of pairs to return, one positive whole number (default
#'   10).
#' @param meth \code{"spearman"} (default) or \code{"pearson"}.
#' @param minabs Report only pairs whose absolute correlation reaches this, a
#'   number in \eqn{[0, 1]} (default 0).
#' @param x A result of \code{topcors()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("topcors", "data.frame")} with
#'   \code{x} and \code{y} (character, column names), \code{correlation}
#'   (numeric) and \code{abscor} (numeric), strongest first, at most \code{n}
#'   rows; the method is in \code{attr(, "method")}.
#' @seealso \code{\link{cortbl}}
#' @examples
#' set.seed(1)
#' n <- 200
#' df <- data.frame(
#'   a = stats::rnorm(n),
#'   b = stats::rnorm(n),
#'   grade = sample(letters[1:3], n, TRUE)
#' )
#' df$c <- df$a * 2 + stats::rnorm(n, sd = 0.1) # strongly related to a
#' df$d <- exp(df$a) # monotone but curved
#'
#' topcors(df)
#'
#' # The curved pair is ranked correctly by Spearman and understated by
#' # Pearson, which is why Spearman is the default.
#' tc <- topcors(df, meth = "spearman")
#' tc[tc$x == "a" & tc$y == "d", "correlation"]
#' tp <- topcors(df, meth = "pearson")
#' tp[tp$x == "a" & tp$y == "d", "correlation"]
#'
#' # Filter to the pairs worth looking at, and keep the top two.
#' topcors(df, n = 2, minabs = 0.5)
#'
#' # a threshold above 1 is an error
#' try(topcors(df, minabs = 2))
#' @export
topcors <- function(data, n = 10L, meth = c("spearman", "pearson"),
                    minabs = 0) {
  f <- "topcors"
  meth <- .ckarg(meth, c("spearman", "pearson"), "meth", f)
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  num <- names(data)[vapply(data, is.numeric, logical(1))]
  if (length(num) < 2L) {
    .ckstop(f, "need at least two numeric columns to correlate")
  }
  n <- .ckpos(n, "n", f)
  minabs <- .ckn(minabs, "minabs", f)
  if (length(minabs) != 1L || is.na(minabs) || minabs < 0 || minabs > 1) {
    .ckstop(f, "`minabs` must be a single value in [0, 1]")
  }
  r <- .Call(C_clab_r_pcor, .nummat(data, num), meth == "pearson", 3L)
  pr <- utils::combn(num, 2L)
  v <- r[[1L]]
  out <- data.frame(x = pr[1L, ], y = pr[2L, ], correlation = v,
                    abscor = abs(v), stringsAsFactors = FALSE)
  out <- out[!is.na(out$correlation) & out$abscor >= minabs, , drop = FALSE]
  out <- utils::head(out[order(-out$abscor), , drop = FALSE], n)
  rownames(out) <- NULL
  attr(out, "method") <- meth
  class(out) <- c("topcors", "data.frame")
  out
}

#' @rdname topcors
#' @export
print.topcors <- function(x, ...) {
  cat("Top correlations (", attr(x, "method"), ")\n", sep = "")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(d, "method") <- NULL
  if (nrow(d)) {
    print(d, row.names = FALSE, digits = 3)
  } else {
    cat("  (none above the threshold)\n")
  }
  invisible(x)
}

#' Drop empty or constant columns and rows
#'
#' \code{drpempty()} removes rows or columns that are entirely missing.
#' \code{drpconst()} removes columns that hold a single distinct value. Use
#' them after reading a spreadsheet export, which often carries blank padding
#' rows, empty columns and columns that repeat one value (a year, a
#' jurisdiction) on every row.
#'
#' @details
#' A constant column carries no information and breaks anything that scales
#' by variance, so it is worth removing -- but it is also a FINDING. A column
#' that was informative in an earlier release and is constant in a fresh one
#' means the source changed, so check that before deleting it and moving on.
#'
#' \code{drpempty()} removes all-\code{NA} columns first, then all-\code{NA}
#' rows among the remaining columns, and resets the row names. An empty
#' string is a value, not missing. \code{drpconst()} counts the distinct
#' non-missing values of each column (with \code{naasval = TRUE}, \code{NA}
#' counts as a value), so an all-\code{NA} column is also dropped by default.
#' Neither function modifies anything in place.
#'
#' Errors: \code{data} not a data frame (or with missing or duplicated
#' names); \code{which} other than \code{"rows"} and/or \code{"cols"}.
#'
#' @param data A data frame whose empty rows or columns are dropped.
#' @param which \code{"rows"}, \code{"cols"}, or both (the default,
#'   \code{c("rows", "cols")}).
#' @param naasval Logical: treat \code{NA} as a distinct value, so a column of
#'   \code{NA}s plus one real value counts as two (default \code{FALSE}).
#' @param quiet Logical: \code{TRUE} skips the summary of what was dropped
#'   (rows and columns before and after, the columns removed), which is
#'   otherwise sent through \code{message()}; default \code{FALSE}.
#' @return The data frame, with the offending rows or columns removed. The
#'   names of the columns that were dropped are attached as the
#'   \code{"dropped"} attribute (a character vector, empty when none);
#'   dropped rows are not listed.
#' @seealso \code{\link{profcols}}, which reports \code{n_distinct} without
#'   removing anything.
#' @examples
#' df <- data.frame(
#'   keep = c(1, 2, NA),
#'   allna = c(NA, NA, NA),
#'   constant = c(7, 7, 7)
#' )
#'
#' drpempty(df)
#' attr(drpempty(df), "dropped")
#'
#' drpconst(df)
#' attr(drpconst(df), "dropped")
#'
#' # Rows only.
#' drpempty(data.frame(a = c(1, NA), b = c(2, NA)), which = "rows")
#'
#' # A column of NAs plus one value is constant by default, and not when
#' # NA is treated as a value of its own.
#' x <- data.frame(v = c(NA, NA, 5))
#' ncol(drpconst(x))
#' ncol(drpconst(x, naasval = TRUE))
#'
#' # only data frames
#' try(drpempty(c(1, NA)))
#' @name drp
#' @export
drpempty <- function(data, which = c("rows", "cols"), quiet = FALSE) {
  f <- "drpempty"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  which <- .ckargs(which, c("rows", "cols"), "which", f)
  quiet <- .ckflag(quiet, "quiet", f)
  before <- dim(data)
  dropped <- character(0)
  if ("cols" %in% which && ncol(data) > 0L) {
    empty <- vapply(data, function(v) all(is.na(v)), logical(1))
    dropped <- names(data)[empty]
    data <- data[, !empty, drop = FALSE]
  }
  if ("rows" %in% which && nrow(data) > 0L && ncol(data) > 0L) {
    na <- vapply(data, is.na, logical(nrow(data)))
    if (is.null(dim(na))) na <- matrix(na, nrow = nrow(data))
    keep <- .Call(C_clab_r_rowany, na) > 0L
    data <- data[keep, , drop = FALSE]
    rownames(data) <- NULL
  }
  attr(data, "dropped") <- dropped
  if (!quiet) {
    .dnote("Empty rows and columns dropped", c(
      sprintf("rows %d -> %d, columns %d -> %d", before[1L], nrow(data),
              before[2L], ncol(data)),
      if (length(dropped)) {
        paste("columns dropped (all missing):", toString(dropped))
      }
    ))
  }
  data
}

#' @rdname drp
#' @export
drpconst <- function(data, naasval = FALSE, quiet = FALSE) {
  f <- "drpconst"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  quiet <- .ckflag(quiet, "quiet", f)
  if (ncol(data) == 0L) {
    attr(data, "dropped") <- character(0)
    return(data)
  }
  const <- vapply(
    data,
    function(v) {
      u <- if (isTRUE(naasval)) unique(v) else unique(v[!is.na(v)])
      length(u) <= 1L
    },
    logical(1)
  )
  dropped <- names(data)[const]
  out <- data[, !const, drop = FALSE]
  attr(out, "dropped") <- dropped
  if (!quiet) {
    .dnote("Constant columns dropped", c(
      sprintf("columns %d -> %d", ncol(data), ncol(out)),
      if (length(dropped)) {
        paste("columns dropped (one value):", toString(dropped))
      }
    ))
  }
  out
}

#' Tidy the column names of an ingested table
#'
#' Converts column names to a consistent, syntactically valid form:
#' transliterated to ASCII, non-alphanumerics collapsed to a single
#' separator, and duplicates disambiguated with a numeric suffix. Use it
#' right after reading an open-data file, so the rest of the script can refer
#' to columns without backticks.
#'
#' @details
#' Open-data extracts arrive with names like \code{"Total Population
#' (2021)"} and \code{"\% change"}, which need backticks everywhere and break
#' silently when a re-release renames \code{"\% change"} to \code{"\%
#' Change"} or adds a trailing space. Normalising once at ingestion makes the
#' schema stable against that.
#'
#' The steps, in order: transliterate to ASCII with [toascii()], which
#' gives the same letters on every system; replace
#' \code{\%} by \code{pct}, \code{&} by \code{and} and \code{#} by \code{n} as
#' words; split camelCase at a lower-to-upper boundary; turn every run of
#' other characters into one word break; an empty name becomes \code{x}; join
#' the words in the chosen case; prefix \code{x} to a name starting with a
#' digit; and give repeated names the suffixes \code{<sep>2}, \code{<sep>3},
#' ... in order of appearance (the first keeps its name).
#'
#' Record the mapping with the data. Cleaning names changes what a
#' downstream script must refer to, so an unrecorded cleaning is itself a
#' reproducibility hazard -- the returned object carries the original names
#' in its \code{"orignms"} attribute for exactly that.
#'
#' @param data A data frame, or a character vector of names.
#' @param case \code{"snake"} (default, \code{total_pop}),
#'   \code{"lower_camel"} (\code{totalPop}), \code{"upper_camel"}
#'   (\code{TotalPop}), \code{"screaming_snake"} (\code{TOTAL_POP}), or
#'   \code{"none"} to normalise separators only and keep the letters' case.
#' @param sep Separator, a string, used between words for \code{"snake"},
#'   \code{"screaming_snake"} and \code{"none"}, and before the number that
#'   disambiguates a duplicate in every case (default \code{"_"}).
#' @param quiet Logical: \code{TRUE} skips the list of the names that
#'   changed, which is otherwise sent through \code{message()} (so
#'   \code{suppressMessages()} also silences it); default \code{FALSE}.
#' @return For a data frame, the data frame with new names and an
#'   \code{"orignms"} attribute (character, the names before cleaning); for a
#'   character vector, the cleaned character vector.
#' @examples
#' clncolnms(c("Total  Population (2021)", "% change", "Ville / City",
#'             "dup", "dup"))
#'
#' # Applied to a data frame, with the original names retained.
#' df <- data.frame(`Total Pop` = 1:2, `% change` = 3:4, check.names = FALSE)
#' cleaned <- clncolnms(df)
#' names(cleaned)
#' attr(cleaned, "orignms")
#'
#' # Other cases.
#' clncolnms(c("Total Pop"), case = "lower_camel")
#' clncolnms(c("Total Pop"), case = "upper_camel")
#' clncolnms(c("Total Pop"), case = "screaming_snake")
#'
#' # a different separator, and a name that starts with a digit
#' clncolnms(c("2021 count", "Rate per 1000"), sep = ".")
#'
#' # an unknown case is an error
#' try(clncolnms("a b", case = "kebab"))
#' @export
clncolnms <- function(
  data,
  case = c("snake", "lower_camel", "upper_camel", "screaming_snake", "none"),
  sep = "_",
  quiet = FALSE
) {
  f <- "clncolnms"
  .cknofn(data, "data", f)
  case <- .ckarg(
    case,
    c("snake", "lower_camel", "upper_camel", "screaming_snake", "none"),
    "case", f
  )
  .cnsep(sep, f)
  quiet <- .ckflag(quiet, "quiet", f)
  .cndata(data, f)
  nms <- if (is.data.frame(data)) names(data) else as.character(data)
  .cktxt(nms, "data", f)
  if (length(nms) == 0L) {
    return(if (is.data.frame(data)) data else character(0))
  }
  out <- .cndup(.cnfix(nms, case, sep), sep)
  if (!quiet) {
    ch <- out != nms
    .dnote(
      "Column names cleaned",
      if (any(ch)) {
        sprintf("%s -> %s", nms[ch], out[ch])
      } else {
        "no name changed"
      }
    )
  }
  if (is.data.frame(data)) {
    names(data) <- out
    attr(data, "orignms") <- nms
    return(data)
  }
  out
}

#' Multivariate outliers by Mahalanobis distance
#'
#' Finds the rows of a numeric table that are unusual as a combination of
#' values, not just on one variable at a time, and gives each row a
#' distance, a p-value and an outlier flag. Use it before modelling to find
#' records that are probably data-entry or linkage errors.
#'
#' @details
#' Ranks rows by how far they sit from the centre of the data ONCE THE
#' CORRELATIONS ARE ACCOUNTED FOR, which is what a per-column check cannot
#' do: a row can be unremarkable on every variable separately and still be
#' impossible jointly -- a person 1.5 m tall weighing 140 kg is inside both
#' marginal ranges and outside the cloud.
#'
#' With centre \eqn{c} and scatter matrix \eqn{S}, the squared distance of
#' row \eqn{x} is \deqn{d^2 = (x - c)^\top S^{-1} (x - c).} Under
#' multivariate normality the squared distance is chi-square on \code{p}
#' degrees of freedom (the number of numeric columns), which gives the
#' p-value and the cut-off.
#'
#' \code{robust = TRUE} centres each column on its median and scales it by
#' its MAD (the 1.4826-scaled median absolute deviation; a column whose MAD
#' is zero falls back to its standard deviation), then takes \eqn{S} as the
#' sample covariance of the standardised data. \code{robust = FALSE} uses the
#' column means and the sample covariance. The robust centre and scales
#' resist outliers; the correlation shape in \eqn{S} is still a classical
#' estimate, so this is not a high-breakdown estimator such as the minimum
#' covariance determinant. The choice matters more than it sounds: outliers
#' inflate the very covariance used to judge them, so with several of them
#' the classical distance hides exactly the rows it is meant to find (the
#' masking effect).
#'
#' Exactly collinear columns are refused rather than repaired: the distance
#' is undefined there, and flooring the covariance's eigenvalues to make it
#' invertible would return numbers governed by the floor rather than by the
#' data. Drop the redundant column first -- \code{\link{topcors}} will
#' identify it.
#'
#' Errors: \code{data} neither a data frame nor a matrix; no numeric columns;
#' no more complete rows than columns; collinear columns (smallest eigenvalue
#' of \eqn{S} at most \eqn{10^{-8}} times the largest); \code{alpha} outside
#' (0, 1).
#'
#' @param data A data frame or numeric matrix; non-numeric columns are
#'   dropped. Rows with any missing value are skipped, and reported as
#'   \code{NA}.
#' @param alpha Significance level for the \code{outlier} flag (default
#'   0.001, deliberately strict: at 0.05 one row in twenty is flagged by
#'   construction).
#' @param robust Use the median/MAD centre and scale (default \code{TRUE}).
#' @param x A result of \code{mahaoutl()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("mahaoutl", "data.frame")} with
#'   \code{row}, \code{distance} (the square root of the squared Mahalanobis
#'   distance), \code{p_value}, \code{logpvalue} (the log of the chi-square
#'   upper tail, which keeps its information where \code{p_value} underflows
#'   to 0) and \code{outlier} (logical), ordered by descending distance, with
#'   skipped rows (\code{NA} distance) last. \code{row} is the row number in
#'   \code{data}. The attributes \code{robust} and \code{alpha} record the
#'   settings. The print method shows the number of flagged rows and the ten
#'   most distant.
#' @references Mahalanobis, P. C. (1936). On the generalised distance in
#'   statistics. Proceedings of the National Institute of Sciences of India
#'   2(1), 49-55.
#' @seealso \code{\link{topcors}}, \code{\link{mcartst}}.
#' @examples
#' set.seed(1)
#' n <- 200
#' df <- data.frame(height = stats::rnorm(n, 170, 10)) # simulated
#' df$weight <- df$height * 0.5 + stats::rnorm(n, 0, 5)
#'
#' # A row that is ordinary on each variable but unlikely jointly:
#' # short and heavy.
#' df[1, ] <- list(height = 155, weight = 100)
#'
#' out <- mahaoutl(df)
#' head(out, 3)
#'
#' # It is flagged, although both values lie inside the other rows' ranges.
#' out$row[1] == 1
#' range(df$height[-1])
#' range(df$weight[-1])
#'
#' # A cluster of bad rows drags the mean toward itself; the median and
#' # MAD move less, so the robust distance of such a row is the larger.
#' bad <- df
#' bad[2:21, ] <- list(height = 150, weight = 110)
#' dclassical <- mahaoutl(bad, robust = FALSE)
#' drobust <- mahaoutl(bad, robust = TRUE)
#' dclassical$distance[dclassical$row == 2] <
#'   drobust$distance[drobust$row == 2]
#'
#' # A duplicated column has no distance defined, and is refused.
#' dup <- df
#' dup$height2 <- dup$height
#' try(mahaoutl(dup))
#' @export
mahaoutl <- function(data, alpha = 0.001, robust = TRUE) {
  f <- "mahaoutl"
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) {
    .ckstop(f, "`data` must be a data frame or a numeric matrix")
  }
  .ckflag(robust, "robust", f)
  num <- vapply(data, is.numeric, logical(1))
  if (!any(num)) .ckstop(f, "`data` has no numeric columns")
  X <- .nummat(data, names(data)[num])
  p <- ncol(X)
  alpha <- .mhalpha(alpha, f)
  ok <- stats::complete.cases(X)
  if (sum(ok) < p + 1L) {
    .ckstop(f, "need more complete rows than columns to estimate a covariance")
  }
  r <- .Call(C_clab_r_maha, X[ok, , drop = FALSE], robust, p)
  if (is.null(r)) {
    .ckstop(
      f,
      paste0(
        "the columns are collinear, so no Mahalanobis distance is ",
        "defined. Drop the redundant column(s) -- drpconst() and ",
        "topcors() will find them"
      )
    )
  }
  out <- data.frame(row = which(ok), distance = r[[1L]], p_value = r[[2L]],
                    logpvalue = r[[3L]], stringsAsFactors = FALSE)
  out$outlier <- out$p_value < alpha
  sk <- which(!ok)
  if (length(sk)) {
    out <- rbind(out, data.frame(row = sk, distance = NA_real_,
                                 p_value = NA_real_, logpvalue = NA_real_,
                                 outlier = NA))
  }
  out <- out[order(-out$distance, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "robust") <- robust
  attr(out, "alpha") <- alpha
  class(out) <- c("mahaoutl", "data.frame")
  out
}

#' @rdname mahaoutl
#' @export
print.mahaoutl <- function(x, ...) {
  cat(sprintf(
    "Mahalanobis outliers (%s): %d of %d rows beyond alpha = %s\n",
    if (isTRUE(attr(x, "robust"))) "robust" else "classical",
    sum(x$outlier, na.rm = TRUE), nrow(x), format(attr(x, "alpha"))
  ))
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(d, "robust") <- attr(d, "alpha") <- NULL
  print(utils::head(d, 10L), row.names = FALSE, digits = 4)
  invisible(x)
}

#' Runs of consecutive missing values
#'
#' Finds the maximal stretches of consecutive \code{NA} in each column, with
#' where each begins and how long it is.
#'
#' Row order carries meaning in a table far more often than people allow for
#' -- a time series, an ordered export, a paginated download -- and a long
#' unbroken run of missingness means something different from the same count
#' scattered about. A run says an instrument was down, a page failed to
#' fetch, or a period was never collected; scattered gaps say individual
#' records failed. The rate cannot distinguish them.
#'
#' @details
#' Each column's missingness indicator \code{is.na(x)} is run-length encoded;
#' every run of \code{TRUE} with length at least \code{minrun} is reported.
#' Rows are taken in their current order, so sort the data by time (or by
#' download order) first. A list column counts an element as missing only
#' when \code{is.na()} says so. Errors: \code{data} not a data frame or
#' without columns, or with duplicated or empty names; \code{minrun} not a
#' positive whole number.
#'
#' @param data A data frame, in the row order that carries meaning.
#' @param minrun Whole number, default 2: report only runs at least this long
#'   (a run of 1 is an isolated gap).
#' @param x A result of \code{missruns()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("missruns", "data.frame")} with
#'   \code{column} (character), \code{start} and \code{end} (integer row
#'   numbers, inclusive) and \code{length} (integer), longest first and then
#'   by column name. Zero rows when there are no qualifying runs.
#' @seealso \code{\link{misspat}} for which columns are missing together,
#'   \code{\link{misssmry}} for the rates.
#' @examples
#' # One long outage and two isolated gaps, with the same total count.
#' df <- data.frame(
#'   outage = c(1, 2, NA, NA, NA, NA, 7, 8),
#'   scattered = c(1, NA, 3, 4, NA, 6, NA, NA)
#' )
#' sum(is.na(df$outage)) == sum(is.na(df$scattered))
#'
#' missruns(df)
#'
#' # Isolated gaps too, by lowering the threshold.
#' missruns(df, minrun = 1)
#'
#' # A complete column has no runs.
#' missruns(data.frame(x = 1:5))
#'
#' try(missruns(df, minrun = 0))
#' @export
missruns <- function(data, minrun = 2L) {
  f <- "missruns"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  if (ncol(data) == 0L) .ckstop(f, "`data` has no columns")
  minrun <- .ckpos(minrun, "minrun", f)
  rows <- lapply(names(data), function(nm) {
    r <- .Call(C_clab_r_nrun, is.na(data[[nm]]), minrun)
    if (!length(r[[1L]])) return(NULL)
    data.frame(column = nm, start = r[[1L]], end = r[[2L]],
               length = r[[2L]] - r[[1L]] + 1L, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(column = character(0), start = integer(0),
                      end = integer(0), length = integer(0),
                      stringsAsFactors = FALSE)
  }
  if (nrow(out)) out <- out[order(-out$length, out$column), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("missruns", "data.frame")
  out
}

#' @rdname missruns
#' @export
print.missruns <- function(x, ...) {
  cat("Runs of consecutive missing values:", nrow(x), "\n")
  if (nrow(x)) print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

#' Text map of where the missing values are
#'
#' Draws the missingness of a whole table as a grid, one character per cell
#' block, a console map that needs no graphics device, so it works over a
#' terminal connection, in a log, and in a plain-text summary.
#'
#' Rows are binned so the map fits \code{height} lines; a block is drawn at
#' the shade its missing proportion falls in. Seeing the table at once is the
#' point: a diagonal band, a block of rows, or a single ragged column are all
#' instantly recognisable shapes that a column of percentages is not.
#'
#' @details
#' The rows are split into at most \code{height} consecutive blocks of nearly
#' equal size; each map line is one block, labelled with its first row
#' number, and each map column is one column of \code{data}. A cell's shade is
#' its block's proportion missing rounded to the nearest quarter: blank (0),
#' light, medium, dark, full (1). Where the console cannot show Unicode the
#' shades are \code{" .-+#"}. Column names are written vertically above the
#' map, cut to \code{width} characters. A table with no rows prints
#' \code{"(no rows)"}. Errors: \code{data} not a data frame or without
#' columns; \code{height} not a positive whole number.
#'
#' @param data A data frame whose missing cells are drawn as a grid.
#' @param height Whole number, default 20: the maximum number of lines (row
#'   blocks) in the map.
#' @param width Whole number, default 12: column names are cut to this many
#'   characters in the vertical legend.
#' @return A character vector of the map's lines, invisibly; the map is
#'   printed with \code{cat()} as a side effect.
#' @seealso \code{\link{misspat}}, \code{\link{missruns}}
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   complete = 1:100,
#'   block = c(rep(NA, 30), 31:100),
#'   scattered = ifelse(stats::runif(100) < 0.3, NA, 1), # simulated gaps
#'   mostlygone = c(1:10, rep(NA, 90))
#' )
#' missgrid(df)
#'
#' # a coarser map with short labels
#' missgrid(df, height = 5, width = 4)
#'
#' try(missgrid(df, height = 0))
#' @export
missgrid <- function(data, height = 20L, width = 12L) {
  f <- "missgrid"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  if (ncol(data) == 0L) .ckstop(f, "`data` has no columns")
  height <- .ckpos(height, "height", f)
  width <- .mgwidth(width, f)
  uni <- .uniok()
  h <- if (uni) "\u2500" else "-"
  rule <- function(t = NULL) {
    if (is.null(t)) return(strrep(h, 66L))
    l <- paste0(strrep(h, 2L), " ", t, " ")
    paste0(l, strrep(h, max(0L, 66L - nchar(l))))
  }
  shades <- .mgshade(uni)
  n <- nrow(data)
  if (n == 0L) {
    lines <- c(rule("Missingness map"), "  (no rows)", rule())
    cat(lines, sep = "\n")
    return(invisible(lines))
  }
  nm <- substring(names(data), 1L, width)
  na <- vapply(data, is.na, logical(n))
  if (is.null(dim(na))) na <- matrix(na, nrow = n)
  g <- .Call(C_clab_r_mgrid, na, min(height, n))
  head <- vapply(seq_len(max(nchar(nm))), function(k) {
    # a name shorter than k keeps its column with a blank
    l <- substring(nm, k, k)
    l[l == ""] <- " "
    paste0("  ", strrep(" ", 6L), " ", paste(l, collapse = ""))
  }, character(1))
  rows <- vapply(seq_along(g[[1L]]), function(i) {
    sprintf("  %6s %s", paste0(g[[1L]][i], ""),
            paste(shades[1L + g[[2L]][i, ]], collapse = ""))
  }, character(1))
  lines <- c(
    rule("Missingness map"), head, rows, "",
    paste0("  legend: '", shades[1L], "' none  '", shades[5L],
           "' all missing"),
    rule()
  )
  cat(lines, sep = "\n")
  cat("\n")
  invisible(lines)
}

#' Little's test for data missing completely at random
#'
#' Tests the MCAR assumption: that whether a value is missing is unrelated to
#' any value in the data, observed or not.
#'
#' The assumption matters because it is what licenses the easy options.
#' Dropping incomplete rows is unbiased under MCAR and biased otherwise; mean
#' imputation understates variance under MCAR and distorts the centre as well
#' otherwise. A small p-value here says the easy options are not available.
#' Run it on the numeric columns of a table before deciding how to handle its
#' missing values.
#'
#' @details
#' Rows are grouped by their pattern of missingness. If missingness is
#' unrelated to the values, then each pattern's observed-variable means should
#' agree with the overall estimates, up to sampling error. The statistic sums
#' those disagreements, \deqn{d^2 = \sum_j n_j (\bar{x}_j - \hat{\mu}_j)'
#' \hat{\Sigma}_j^{-1} (\bar{x}_j - \hat{\mu}_j),} over the variables
#' observed in pattern \eqn{j}, and is compared with a chi-square
#' distribution on \eqn{\sum_j p_j - p} degrees of freedom.
#'
#' The overall \eqn{\hat\mu} and \eqn{\hat\Sigma} are the MAXIMUM-LIKELIHOOD
#' estimates under multivariate normality WITH the missing data, obtained by
#' expectation-maximisation -- not the complete-case estimates, which would
#' already embed the bias the test is looking for.
#'
#' A large p-value is NOT evidence that the data are MCAR; it is a failure to
#' detect a departure, and the test has little power on small samples or
#' with many patterns. The test also assumes multivariate normality, so on
#' markedly non-normal columns a rejection may be telling you about the
#' distribution rather than the missingness. Neither this test nor any other
#' can distinguish missing-at-random from missing-NOT-at-random, because that
#' distinction depends on the values that were never observed.
#'
#' Collinear or constant columns are refused rather than worked around: the
#' likelihood is degenerate there, and the EM step's eigenvalue floor would
#' otherwise return a statistic governed by that floor instead of by the
#' data. Columns that are entirely missing and rows with no observed value are
#' dropped before the test. If EM has not converged after \code{maxit}
#' iterations the last estimates are used; \code{iterations} then equals
#' \code{maxit}.
#'
#' Errors: \code{data} neither a data frame nor a matrix; fewer than two
#' numeric columns, or fewer than two usable rows or columns after the drops;
#' collinear or constant columns (smallest covariance eigenvalue at most
#' \eqn{10^{-8}} times the largest); EM estimates that cannot be inverted.
#' Non-numeric columns are dropped with a warning.
#'
#' @param data A data frame or numeric matrix. Non-numeric columns are
#'   dropped with a warning, since the statistic is defined on moments.
#' @param maxit Whole number, default 500: the maximum number of EM
#'   iterations.
#' @param tol Number, default 1e-7: EM stops when the largest absolute change
#'   in any mean or covariance entry between iterations falls below it.
#' @param x A result of \code{mcartst()}.
#' @inheritParams freqtbl
#' @return A list of class \code{c("mcartst", "list")}: \code{statistic},
#'   \code{df}, \code{logpvalue} (the log of \code{p_value}, finite where
#'   \code{p_value} underflows to 0), \code{p_value}, \code{npat} (distinct
#'   missingness patterns), \code{nused} (rows in the patterns that entered
#'   the statistic), \code{nvars} (columns used), \code{iterations} (EM
#'   iterations), \code{mu} (named ML means), \code{sigma} (ML covariance
#'   matrix), and \code{note} (a caveat when the test is degenerate or a
#'   pattern was skipped, else \code{NA}). With complete data \code{df} is 0
#'   and \code{p_value} is \code{NA} -- there is nothing to test.
#' @references Little, R. J. A. (1988). A test of missing completely at
#'   random for multivariate data with missing values. Journal of the
#'   American Statistical Association 83(404), 1198-1202.
#'   \doi{10.1080/01621459.1988.10478722}
#'
#'   Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum likelihood
#'   from incomplete data via the EM algorithm. Journal of the Royal
#'   Statistical Society B 39(1), 1-38.
#' @seealso \code{\link{misspat}} for the patterns themselves,
#'   \code{\link{misssmry}} for the rates.
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- stats::rnorm(n) # simulated
#' y <- x + stats::rnorm(n)
#'
#' # Missing completely at random: a coin flip decides, so the test
#' # should not reject.
#' mcar <- data.frame(x = x, y = y)
#' mcar$y[sample(n, 90)] <- NA
#' mcartst(mcar)
#'
#' # Missing depending on the OTHER, observed variable: not MCAR, and
#' # detectable, because the pattern's mean of x is shifted.
#' mar <- data.frame(x = x, y = y)
#' mar$y[x > 0.4] <- NA
#' mcartst(mar)
#'
#' # Complete data has one pattern and nothing to test.
#' mcartst(data.frame(a = x, b = y))$df
#'
#' # A duplicated column makes the likelihood degenerate, and is refused.
#' dup <- mcar
#' dup$x2 <- dup$x
#' try(mcartst(dup))
#' @export
mcartst <- function(data, maxit = 500L, tol = 1e-7) {
  f <- "mcartst"
  ctl <- .mcctl(maxit, tol, f)
  X <- .mcmat(data, f)
  r <- .Call(C_clab_r_mcar, X, as.integer(ctl$maxit), as.numeric(ctl$tol))
  if (is.null(r)) {
    .ckstop(
      f,
      paste0(
        "the columns are collinear (or one is constant), so the ",
        "multivariate normal likelihood Little's test is built on is ",
        "degenerate. Drop the redundant column(s) first -- ",
        "topcors() and drpconst() will find them"
      )
    )
  }
  o <- r[[1L]]
  df <- o[2L]
  structure(
    list(
      statistic = o[1L],
      df = df,
      logpvalue = if (df > 0) o[4L] else NA_real_,
      p_value = if (df > 0) o[3L] else NA_real_,
      npat = as.integer(o[5L]),
      nused = as.integer(o[6L]),
      nvars = ncol(X),
      iterations = as.integer(o[7L]),
      mu = stats::setNames(r[[2L]], colnames(X)),
      sigma = r[[3L]],
      note = .mcnote(o)
    ),
    class = c("mcartst", "list")
  )
}

#' @rdname mcartst
#' @export
print.mcartst <- function(x, ...) {
  v <- if (is.na(x$p_value)) {
    "no test possible"
  } else if (x$p_value < 0.05) {
    "MCAR rejected: the missingness is related to the data"
  } else {
    "no departure from MCAR detected"
  }
  cat(
    "Little's MCAR test: ", v, "\n",
    "  Statistic ", formatC(x$statistic, format = "f", digits = 3),
    " on ", x$df, " df, p = ", format(x$p_value, digits = 3), "\n",
    "  Patterns ", x$npat, ", variables ", x$nvars, ", rows used ", x$nused,
    ", EM iterations ", x$iterations, "\n",
    sep = ""
  )
  if (!is.na(x$note)) cat("  Note: ", x$note, "\n", sep = "")
  invisible(x)
}

#' Which columns are missing together
#'
#' Counts the distinct PATTERNS of missingness across rows -- which
#' combination of columns is missing in each row -- rather than the
#' per-column rates \code{\link{profcols}} reports. Use it before deciding how
#' to handle gaps (drop rows, impute, model separately), because the right
#' choice depends on whether columns are missing together.
#'
#' @details
#' The distinction decides what to do about the gaps. Two columns each 30\%
#' missing at random need different handling from two columns 30\% missing in
#' THE SAME rows -- the second is one structural gap (a join that failed, a
#' form section nobody filled in) and often means those rows should be
#' dropped or modelled separately, while the first does not. A per-column
#' rate cannot tell the two apart; this can.
#'
#' Each row is coded as a string with one character per column, \code{X}
#' where \code{is.na()} is true and \code{.} otherwise, and the codes are
#' tabulated. Patterns are sorted by descending frequency and the first
#' \code{maxpat} kept, so \code{pctrows} may not sum to 100 when patterns are
#' cut. A frame with no rows gives an empty table.
#'
#' Errors: \code{data} not a data frame, without columns, or with missing or
#' duplicated names; \code{maxpat} not one positive whole number.
#'
#' @param data A data frame whose missingness patterns are counted.
#' @param maxpat Maximum patterns to return, most frequent first, one
#'   positive whole number (default 20).
#' @param x A result of \code{misspat()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("misspat", "data.frame")}, one row
#'   per pattern: \code{pattern} (character, \code{.} for present and
#'   \code{X} for missing, in column order), \code{n_rows} (integer),
#'   \code{pctrows} (numeric, percent of all rows), \code{n_missing}
#'   (integer, columns missing in that pattern) and \code{columns}
#'   (character, their names comma-separated, \code{""} for the complete
#'   pattern). Carries the column order as the \code{"columns"} attribute.
#' @seealso \code{\link{profcols}} for per-column rates.
#' @examples
#' # Two columns missing in the SAME rows: one structural gap.
#' structural <- data.frame(
#'   id = 1:10,
#'   a = c(rep(NA, 3), 4:10),
#'   b = c(rep(NA, 3), 4:10)
#' )
#' misspat(structural)
#'
#' # The same per-column rates, but missing independently.
#' scattered <- data.frame(
#'   id = 1:10,
#'   a = c(rep(NA, 3), 4:10),
#'   b = c(1:7, rep(NA, 3))
#' )
#' misspat(scattered)
#'
#' # A complete frame has exactly one pattern.
#' misspat(data.frame(x = 1:3, y = 4:6))
#'
#' # keep only the most common pattern
#' misspat(scattered, maxpat = 1)
#'
#' # a frame without columns has no patterns
#' try(misspat(data.frame()))
#' @export
misspat <- function(data, maxpat = 20L) {
  f <- "misspat"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  if (ncol(data) == 0L) .ckstop(f, "`data` has no columns")
  maxpat <- .ckpos(maxpat, "maxpat", f)
  n <- nrow(data)
  out <- data.frame(pattern = character(0), n_rows = integer(0),
                    pctrows = numeric(0), n_missing = integer(0),
                    columns = character(0), stringsAsFactors = FALSE)
  if (n > 0L) {
    na <- vapply(data, is.na, logical(n))
    if (is.null(dim(na))) na <- matrix(na, nrow = n)
    codes <- apply(na, 1L, function(r) {
      paste(ifelse(r, "X", "."), collapse = "")
    })
    lv <- sort(unique(codes))
    r <- .Call(C_clab_r_freq, match(codes, lv), length(lv), 0L,
               as.numeric(n), TRUE, maxpat, FALSE)
    keep <- lv[r[[1L]]]
    flags <- strsplit(keep, "", fixed = TRUE)
    out <- data.frame(
      pattern = keep,
      n_rows = r[[2L]],
      pctrows = r[[3L]],
      n_missing = vapply(flags, function(z) sum(z == "X"), integer(1)),
      columns = vapply(flags, function(z) {
        if (any(z == "X")) toString(names(data)[z == "X"]) else ""
      }, ""),
      stringsAsFactors = FALSE
    )
  }
  attr(out, "columns") <- names(data)
  class(out) <- c("misspat", "data.frame")
  out
}

#' @rdname misspat
#' @export
print.misspat <- function(x, ...) {
  cat("Missingness patterns, columns in pattern order:",
      toString(attr(x, "columns")), "\n")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(d, "columns") <- NULL
  print(d, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Duplicated rows, with their groups
#'
#' Returns the rows that share a combination of \code{cols} with at least one
#' other row, grouped so the duplicates sit together. A duplicate is usually
#' a join that fanned out or a re-release appended rather than replaced,
#' both of which are visible only once the offending rows are in front of
#' you.
#'
#' @details
#' Two rows are duplicates when their values, converted to text, are equal in
#' every column of \code{cols}; two \code{NA}s are equal. Columns are taken
#' by position when \code{cols = NULL}, so a table with repeated column names
#' can still be inspected. Errors when \code{data} is not a data frame, a
#' named column does not exist, or a named column's name is shared by two
#' columns.
#'
#' @param data A data frame whose rows are checked for duplicates.
#' @param cols Character vector of the columns defining a duplicate, or
#'   \code{NULL} (the default) for all of them.
#' @param quiet Logical: \code{TRUE} skips the count of duplicate rows and
#'   groups, which is otherwise sent through \code{message()}; default
#'   \code{FALSE}.
#' @return A data frame: the duplicated rows (every occurrence, including the
#'   first), ordered so each group sits together, with an integer
#'   \code{dupe_count} column giving each group's size and row names reset.
#'   Zero rows (with the same columns) when there are none.
#' @seealso \code{\link{misssmry}}
#' @examples
#' df <- data.frame(
#'   id = c(1, 2, 2, 3, 3, 3),
#'   value = c("a", "b", "b", "c", "d", "c")
#' )
#'
#' # Duplicated on every column.
#' duprows(df)
#'
#' # Duplicated on the identifier alone, which catches more.
#' duprows(df, "id")
#'
#' # No duplicates gives zero rows, not an error.
#' duprows(data.frame(x = 1:3))
#'
#' # a column that does not exist is an error
#' try(duprows(df, "code"))
#' @export
duprows <- function(data, cols = NULL, quiet = FALSE) {
  f <- "duprows"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  quiet <- .ckflag(quiet, "quiet", f)
  if (nrow(data) == 0L) {
    out <- cbind(data, dupe_count = integer(0))
  } else {
    pos <- if (is.null(cols)) {
      seq_along(data)
    } else {
      .cknofn(cols, "cols", f)
      cc <- as.character(cols)
      miss <- setdiff(cc, names(data))
      if (length(miss)) .ckstop(f, "no such column(s): %s", toString(miss))
      tw <- cc[cc %in% names(data)[duplicated(names(data))]]
      if (length(tw)) {
        .ckstop(
          f,
          paste0(
            "`data` has two columns named %s: give the columns by ",
            "position, or rename one"
          ),
          tw[[1L]]
        )
      }
      match(cc, names(data))
    }
    key <- do.call(paste, c(lapply(pos, function(i) {
      as.character(data[[i]])
    }), sep = "\r"))
    lv <- sort(unique(key))
    r <- .Call(C_clab_r_dups, match(key, lv), length(lv))
    out <- data[r[[1L]], , drop = FALSE]
    out$dupe_count <- r[[2L]]
    rownames(out) <- NULL
  }
  if (!quiet) {
    .dnote("Duplicated rows", sprintf("%d row(s) of %d", nrow(out),
                                      nrow(data)))
  }
  out
}

#' Missingness in one line per question
#'
#' The scalar summaries of missingness: how much of the table is missing, how
#' many rows are complete, and how many columns have any or only missing
#' values. Use it as a first look at a new extract, or to report data
#' completeness in a methods section; \code{\link{misspat}} then shows which
#' columns are missing together.
#'
#' @details
#' A cell is missing when \code{is.na()} is \code{TRUE} (so \code{NaN}
#' counts, and an empty string does not). A row is complete when no cell in
#' it is missing (\code{stats::complete.cases()}). Percentages are of all
#' cells and of all rows. A table with no rows gives zeros. Columns are taken
#' by position, so duplicate or empty column names are fine. Errors when
#' \code{data} is not a data frame or has no columns.
#'
#' @param data A data frame with at least one column.
#' @return A named numeric vector: \code{n_rows}, \code{n_cols},
#'   \code{n_missing} (missing cells), \code{pct_missing} (percent of cells),
#'   \code{n_complete_rows}, \code{pct_complete_rows} (percent of rows),
#'   \code{n_cols_any_missing} (columns with at least one missing value) and
#'   \code{n_cols_all_missing} (columns with no value at all).
#' @seealso \code{\link{misspat}} for which columns are missing together,
#'   \code{\link{profcols}} for per-column rates.
#' @examples
#' df <- data.frame(a = c(1, NA, 3), b = c(NA, NA, 3), c = 1:3)
#' misssmry(df)
#'
#' # A complete table is all zeros but for its dimensions.
#' misssmry(data.frame(x = 1:3, y = 4:6))
#'
#' # a wholly empty column
#' ms <- misssmry(data.frame(x = 1:4, note = NA))
#' ms[c("pct_missing", "n_cols_all_missing")]
#'
#' # a vector is not a table
#' try(misssmry(c(1, NA, 3)))
#' @export
misssmry <- function(data) {
  f <- "misssmry"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  if (ncol(data) == 0L) .ckstop(f, "`data` has no columns")
  n <- nrow(data)
  na <- vapply(data, function(v) as.vector(is.na(v)), logical(n))
  if (is.null(dim(na))) na <- matrix(na, nrow = n)
  stats::setNames(
    .Call(C_clab_r_msum, na),
    c("n_rows", "n_cols", "n_missing", "pct_missing", "n_complete_rows",
      "pct_complete_rows", "n_cols_any_missing", "n_cols_all_missing")
  )
}

#' Column report for a data frame
#'
#' One row per column: type, missingness, and the summaries that suit the
#' column's type -- mean and standard deviation alongside their robust
#' counterparts (median, MAD) for a numeric column, and the number of
#' distinct levels plus the most common one for a categorical column.
#'
#' Reporting the classical and robust centres side by side is the point:
#' where they disagree, the column has outliers or a heavy tail, and the mean
#' is not describing it. That is visible in one glance here and in no single
#' number.
#'
#' @details
#' A column is numeric when \code{is.numeric()} is \code{TRUE}; every other
#' column (character, factor, logical, \code{Date}, ...) is treated as
#' categorical, and \code{top} is its most frequent non-missing value as text
#' (ties broken by the order of \code{sort(table())}). For a numeric column,
#' over its non-missing values: \code{sd} is the sample standard deviation
#' (\eqn{n - 1}), \code{mad} the median absolute deviation scaled by 1.4826,
#' \code{skewness} the population moment ratio \eqn{m_3 / m_2^{3/2}},
#' \code{n_outliers} the count outside the Tukey fences (1.5 interquartile
#' ranges beyond the quartiles), and the quantiles are type 7. Infinite values
#' are counted in \code{ninf} (and in \code{nneg}, \code{n_distinct}) and left
#' out of every other numeric statistic, which are taken over the finite
#' values only: a column with an infinite value keeps a finite mean and
#' range, and a column with no finite value has \code{NA} statistics. Errors
#' when \code{data} is not a data frame, has no columns, or has duplicate or
#' empty column names.
#'
#' @param data A data frame whose columns are profiled.
#' @param qnts Numeric vector of quantile probabilities to include for
#'   numeric columns, one column each named \code{q} and the percentage
#'   (\code{q25}, \code{q50}, ...); default the quartiles. \code{NULL} or an
#'   empty vector adds none.
#' @param x A result of \code{profcols()}.
#' @inheritParams freqtbl
#' @return A data frame of class \code{c("profcols", "data.frame")}, one row
#'   per column, with the columns \code{column}, \code{type} (the first
#'   class), \code{n}, \code{n_missing}, \code{pct_missing},
#'   \code{n_distinct}, \code{nzero}, \code{nneg}, \code{ninf}, \code{mean},
#'   \code{sd}, \code{median}, \code{mad}, \code{min}, \code{max},
#'   \code{skewness}, \code{n_outliers}, \code{hist} (an \code{\link{inlhist}}
#'   sketch), \code{top} and one per quantile; a statistic that does not apply
#'   to the column's type is \code{NA}. The print method shows the histogram
#'   last.
#' @seealso \code{\link{inlhist}}, \code{\link{freqtbl}}.
#' @examples
#' set.seed(9)
#' df <- data.frame( # simulated columns
#'   clean = stats::rnorm(100),
#'   skewed = c(stats::rnorm(99), 500),
#'   grade = sample(c("a", "b", "c"), 100, TRUE),
#'   gappy = c(rep(NA, 10), stats::runif(90))
#' )
#'
#' profcols(df)
#'
#' # The mean and the median agree on `clean` and disagree sharply on
#' # `skewed`, which is the outlier announcing itself.
#' p <- profcols(df)
#' p[p$column %in% c("clean", "skewed"), c("column", "mean", "median")]
#'
#' # The histogram column shows shape no summary number carries.
#' p[, c("column", "hist")]
#'
#' # Zeros, negatives and infinities are counted separately, because
#' # each breaks a different downstream computation (a log, a square
#' # root, an average).
#' p[, c("column", "nzero", "nneg", "ninf")]
#'
#' # deciles instead of quartiles
#' names(profcols(df, qnts = c(0.1, 0.9)))
#'
#' # a table with no columns has nothing to profile
#' try(profcols(df[, 0]))
#' @export
profcols <- function(data, qnts = c(0.25, 0.5, 0.75)) {
  f <- "profcols"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  .dfnm(data, f)
  if (ncol(data) == 0L) .ckstop(f, "`data` has no columns to profile")
  qnts <- .ckn(qnts, "qnts", f)
  if (anyNA(qnts) || any(qnts < 0 | qnts > 1)) {
    .ckstop(f, "`qnts` must be probabilities in [0, 1]")
  }
  rows <- lapply(names(data), function(nm) .profrow(data[[nm]], nm))
  out <- do.call(rbind, rows)
  if (length(qnts)) {
    qm <- do.call(rbind, lapply(data, function(v) {
      if (is.numeric(v) && any(is.finite(v))) {
        .Call(C_clab_r_q7, as.numeric(v[is.finite(v)]), as.numeric(qnts))
      } else {
        rep(NA_real_, length(qnts))
      }
    }))
    colnames(qm) <- paste0("q", format(100 * qnts, trim = TRUE))
    rownames(qm) <- NULL
    out <- cbind(out, qm)
  }
  rownames(out) <- NULL
  class(out) <- c("profcols", "data.frame")
  out
}

#' @rdname profcols
#' @export
print.profcols <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if ("hist" %in% names(d)) d <- d[c(setdiff(names(d), "hist"), "hist")]
  cat("Column profile:", nrow(d), "columns\n")
  print(d, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Do distinct groups differ? A heterogeneity test before pooling
#'
#' Tests whether the estimates of one indicator for distinct groups (First
#' Nations, Metis and Inuit; Aboriginal and Torres Strait Islander people;
#' any groups with published intervals) differ more than their intervals
#' allow, and says how much of their spread is real: the evidence for
#' reporting each group distinctly rather than as one pooled figure.
#'
#' @details With estimates \eqn{y_g}, standard errors \eqn{s_g} from the
#' published intervals and weights \eqn{w_g = 1/s_g^2}, Cochran's \eqn{Q =
#' \sum w_g (y_g - \bar y_w)^2} has a chi-squared distribution on \eqn{G - 1}
#' degrees of freedom when the groups share one value, and \eqn{I^2 = \max(0,
#' (Q - (G-1))/Q)} is the share of the spread that is real rather than
#' sampling error (Higgins and Thompson 2002). A pooled figure is a weighted
#' mean of different populations whenever \eqn{I^2} is large; a
#' population-weighted total also depends on how many people are in each
#' group, so it is mostly the largest group's figure. Rows that pool groups
#' are left out (\code{pooled}, or labels starting with "Total", "Indigenous
#' identity", "Aboriginal identity", "Total Indigenous" or "Indigenous
#' population"); every other row is included, so drop a reference group
#' yourself, and never test overlapping groups (\code{"alone"} beside
#' \dQuote{alone or in combination}) together: they are not independent. Rows
#' with a quality flag \code{"F"}, \code{"x"}, \code{".."} or \code{"..."}
#' (suppressed) are left out; other flags than those and \code{"E"} are
#' treated as usable, with a warning. An estimate outside its own interval is
#' not compared, with a warning. With two or three groups the test has little
#' power and \eqn{I^2} is imprecise; a large \eqn{p} is not evidence that the
#' groups are the same.
#'
#' @param data A data frame, one row per group.
#' @param grp,est,lower,upper Column names.
#' @param qual Optional column of quality flags; suppressed rows are left out.
#' @param conflvl Confidence level of the published intervals, and the level
#'   of the verdict printed (default 0.95).
#' @param pooled Optional labels of rows that pool other groups.
#' @param x A result of \code{dstntst()}.
#' @inheritParams freqtbl
#' @return A list of class \code{dstntst}: \code{groups} (the estimates used,
#'   with their standard errors), \code{pooled} (the inverse-variance mean)
#'   and \code{poose}, \code{q}, \code{df}, \code{p_value}, \code{i2},
#'   \code{leftout} (the rows not used) and \code{conflvl}. \code{print()}
#'   shows it.
#' @references Higgins, J. P. T. and Thompson, S. G. (2002). Quantifying
#'   heterogeneity in a meta-analysis. Statistics in Medicine 21(11),
#'   1539-1558. \doi{10.1002/sim.1186}
#' @examples
#' # your own table: made-up numbers, with the pooled row named
#' own <- data.frame(
#'   group = c("Aboriginal", "Torres Strait Islander", "Both"),
#'   est = c(24.0, 31.0, 24.6),
#'   lo = c(22.1, 26.2, 22.8),
#'   hi = c(25.9, 35.8, 26.4)
#' )
#' dstntst(own, "group", "est", "lo", "hi", pooled = "Both")
#'
#' # simulated victimisation rates of three police divisions
#' div <- data.frame(
#'   division = c("North", "Central", "South"),
#'   rate = c(4.1, 5.3, 4.4), lo = c(3.6, 4.7, 3.9), hi = c(4.6, 5.9, 4.9)
#' )
#' dstntst(div, "division", "rate", "lo", "hi")
#' @export
dstntst <- function(
  data,
  grp,
  est,
  lower,
  upper,
  qual = NULL,
  conflvl = 0.95,
  pooled = NULL
) {
  f <- "dstntst"
  if (!is.data.frame(data)) .ckstop(f, "`data` must be a data frame")
  grp <- .dtcol(grp, "grp", data, f)
  est <- .dtcol(est, "est", data, f)
  lower <- .dtcol(lower, "lower", data, f)
  upper <- .dtcol(upper, "upper", data, f)
  if (!is.null(qual)) qual <- .dtcol(qual, "qual", data, f)
  conflvl <- .dfcl(conflvl, f)
  z <- stats::qnorm(1 - (1 - conflvl) / 2)
  g <- enc2utf8(trimws(as.character(data[[grp]])))
  .dtpool(pooled, g, f)
  y <- .ckn(data[[est]], "est", f)
  lo <- .ckn(data[[lower]], "lower", f)
  hi <- .ckn(data[[upper]], "upper", f)
  y <- .dtest(y, lo, hi, g, f)
  s <- (hi - lo) / (2 * z)
  q <- .dtqual(if (is.null(qual)) rep("ok", length(g)) else data[[qual]], f)
  pan <- .dtpan(g, pooled)
  use <- !pan & !q %in% c("F", "x", "..", "...") & is.finite(y) &
    is.finite(s) & s > 0
  if (sum(use) < 2L) {
    .ckstop(
      f,
      "the test needs two groups at least with an estimate and an interval"
    )
  }
  if (anyDuplicated(g[use])) .ckstop(f, "a group appears twice")
  r <- .Call(C_clab_r_cochq, as.numeric(y[use]), as.numeric(s[use]))
  structure(
    list(
      groups = data.frame(group = g[use], estimate = y[use], se = s[use],
                          stringsAsFactors = FALSE),
      pooled = r[1L],
      poose = r[2L],
      q = r[3L],
      df = sum(use) - 1L,
      p_value = r[5L],
      i2 = r[6L],
      leftout = g[!use],
      conflvl = conflvl
    ),
    class = "dstntst"
  )
}

#' @rdname dstntst
#' @export
print.dstntst <- function(x, ...) {
  cat(
    "Distinctions test: ", nrow(x$groups), " groups\n",
    sprintf("  Cochran's Q %.2f on %d df, p = %s\n", x$q, x$df,
            format.pval(x$p_value, digits = 3)),
    sprintf("  I-squared: %.0f%% of the variation is beyond chance\n",
            100 * x$i2),
    sprintf("  Pooled: %.2f (se %.3g, inverse-variance mean)\n", x$pooled,
            x$poose),
    if (length(x$leftout)) {
      paste0("  Left out: ", paste(x$leftout, collapse = "; "), "\n")
    },
    sep = ""
  )
  print(x$groups, row.names = FALSE)
  invisible(x)
}

# .rmbl_df_names(): a data frame whose columns all have distinct names.
.dfnm <- function(data, fn, what = "data") {
  if (!is.data.frame(data)) {
    .ckstop(
      fn, "`%s` must be a data frame, not %s", what,
      if (is.null(data)) "NULL" else class(data)[[1L]]
    )
  }
  nm <- names(data)
  if (length(data) && (is.null(nm) || anyNA(nm) || !all(nzchar(nm)))) {
    .ckstop(
      fn,
      "every column of `%s` must have a name (names(%s) has an empty one)",
      what, what
    )
  }
  if (anyDuplicated(nm)) {
    .ckstop(fn, "`%s` has two columns named %s: rename one", what,
            nm[duplicated(nm)][[1L]])
  }
  invisible(TRUE)
}

# .rmbl_num(): numbers, text that reads as numbers, or missing values; a
# whole number when integer.
.ckn <- function(x, what, fn, integer = FALSE) {
  if (!is.null(x) && (!is.atomic(x) || is.factor(x) || is.complex(x) ||
                        is.raw(x))) {
    .ckstop(
      fn, "`%s` must be a%s number, not %s", what,
      if (integer) " whole" else "",
      if (is.data.frame(x)) {
        "a data frame"
      } else if (is.factor(x)) {
        "a factor"
      } else if (is.function(x)) {
        "a function"
      } else {
        class(x)[[1L]]
      }
    )
  }
  if (is.character(x)) {
    .cktxt(x, what, fn)
    v <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & is.na(v) & nzchar(trimws(x)) & x != "NA"
    if (any(bad)) {
      .ckstop(fn, "`%s` must be numeric, not text such as \"%s\"", what,
              x[bad][1L])
    }
    x <- v
  }
  if (is.logical(x) && !all(is.na(x))) {
    .ckstop(fn, "`%s` must be a%s number, not %s", what,
            if (integer) " whole" else "", format(x[!is.na(x)][1L]))
  }
  x <- as.numeric(x)
  if (!integer) {
    return(x)
  }
  frac <- (is.finite(x) & x != floor(x)) | is.infinite(x)
  if (any(frac)) {
    .ckstop(fn, "`%s` must be a whole number, not %s", what,
            format(x[frac][1L]))
  }
  big <- !is.na(x) & abs(x) > .Machine$integer.max
  if (any(big)) {
    .ckstop(fn, "`%s` is too large: %s (the limit is 2147483647)", what,
            format(x[big][1L]))
  }
  as.integer(x)
}

# One positive whole number.
.ckpos <- function(x, what, fn) {
  v <- .ckn(x, what, fn, integer = TRUE)
  if (length(v) != 1L || is.na(v) || v < 1L) {
    .ckstop(fn, "`%s` must be a single positive integer", what)
  }
  v
}

# Data, not a function or an environment.
.cknofn <- function(x, what, fn) {
  if (is.function(x)) .ckstop(fn, "`%s` must be data, not a function", what)
  if (is.environment(x)) {
    .ckstop(fn, "`%s` must be data, not an environment", what)
  }
  invisible(TRUE)
}

# One or more of the choices, matched partially.
.ckargs <- function(arg, choices, what, fn) {
  if (is.null(arg)) {
    return(choices[[1L]])
  }
  i <- if (is.character(arg) && length(arg) && !anyNA(arg)) {
    pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  } else {
    0L
  }
  if (any(i == 0L)) {
    .ckstop(fn, "`%s` must be one or more of %s, not %s", what,
            toString(dQuote(choices, FALSE)),
            paste(utils::head(deparse(arg), 1L), collapse = ""))
  }
  unique(choices[i])
}

# Unicode on this console, unless CRIMLAB_ASCII_ONLY is set.
.uniok <- function() {
  identical(Sys.getenv("CRIMLAB_ASCII_ONLY"), "") &&
    isTRUE(l10n_info()[["UTF-8"]])
}

# A short boxed note through message().
.dnote <- function(title, lines) {
  message(paste(c(paste0("-- ", title), paste0("   ", lines)), collapse = "\n"))
}

# The sparkline of a numeric vector of finite values.
.spark <- function(x, bins) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(strrep(" ", bins))
  }
  uni <- .uniok()
  if (min(x) == max(x)) {
    return(paste0(
      strrep(" ", bins %/% 2L), if (uni) "\u2588" else "#",
      strrep(" ", bins - bins %/% 2L - 1L)
    ))
  }
  cnt <- .Call(C_clab_r_hbin, as.numeric(x), as.integer(bins))
  lev <- if (uni) {
    c("\u2581", "\u2582", "\u2583", "\u2584", "\u2585", "\u2586", "\u2587",
      "\u2588")
  } else {
    c(".", ":", "-", "=", "#")
  }
  idx <- ifelse(cnt == 0L, 0L, pmax(1L, ceiling(cnt / max(cnt) * length(lev))))
  paste(ifelse(idx == 0L, " ", lev[pmax(1L, idx)]), collapse = "")
}

# The numeric columns of a data frame as a double matrix, NA as NaN.
.nummat <- function(data, cols) {
  m <- matrix(NA_real_, nrow(data), length(cols))
  for (j in seq_along(cols)) m[, j] <- as.numeric(data[[cols[j]]])
  m
}

.cnsep <- function(sep, f) {
  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    .ckstop(f, "`sep` must be a single string")
  }
  .cktxt(sep, "sep", f)
}

.cndata <- function(data, f) {
  if (!is.data.frame(data) &&
        (!(is.character(data) || is.factor(data)) || anyNA(data))) {
    .ckstop(f, "`data` must be a data frame or column names, without NA")
  }
}

# clncolnms(): the names in ASCII, split into words and joined in `case`
.cnfix <- function(nms, case, sep) {
  # not iconv(, "ASCII//TRANSLIT"): its letters depend on the system's
  # iconv ("ete" on one, "'et'e" on another for the same accented word)
  out <- toascii(nms)
  out[is.na(out)] <- nms[is.na(out)]
  out <- gsub("%", " pct ", out, fixed = TRUE)
  out <- gsub("&", " and ", out, fixed = TRUE)
  out <- gsub("#", " n ", out, fixed = TRUE)
  out <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", out)
  out <- trimws(gsub("[^A-Za-z0-9]+", " ", out))
  out[!nzchar(out)] <- "x"
  cap <- function(w) {
    paste0(toupper(substring(w, 1L, 1L)), tolower(substring(w, 2L)))
  }
  out <- vapply(strsplit(out, " ", fixed = TRUE), function(w) {
    w <- w[nzchar(w)]
    switch(
      case,
      snake = paste(tolower(w), collapse = sep),
      screaming_snake = paste(toupper(w), collapse = sep),
      lower_camel = paste0(tolower(w[1L]), paste(cap(w[-1L]), collapse = "")),
      upper_camel = paste(cap(w), collapse = ""),
      none = paste(w, collapse = sep)
    )
  }, character(1))
  ifelse(grepl("^[0-9]", out), paste0("x", out), out)
}

# clncolnms(): a repeated name numbered from its second time, sep 2, 3, ...
.cndup <- function(out, sep) {
  dup <- duplicated(out)
  if (any(dup)) {
    for (nm in unique(out[dup])) {
      i <- which(out == nm)
      out[i] <- paste0(nm, c("", paste0(sep, seq_len(length(i) - 1L) + 1L)))
    }
  }
  out
}

.mhalpha <- function(alpha, f) {
  alpha <- .ckn(alpha, "alpha", f)
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    .ckstop(f, "`alpha` must be a single value strictly inside (0, 1)")
  }
  alpha
}

.mgwidth <- function(width, f) {
  width <- .ckn(width, "width", f, integer = TRUE)
  if (length(width) != 1L || is.na(width) || width < 1 || width > 10000) {
    .ckstop(f, "`width` must be one whole number from 1 to 10000")
  }
  width
}

# missgrid(): the five shades from none to all missing
.mgshade <- function(uni) {
  if (uni) {
    c(" ", "\u2591", "\u2592", "\u2593", "\u2588")
  } else {
    c(" ", ".", "-", "+", "#")
  }
}

.mcctl <- function(maxit, tol, f) {
  maxit <- .ckn(maxit, "maxit", f, integer = TRUE)
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1 || maxit > 1e6) {
    .ckstop(f, "`maxit` must be one whole number from 1 to 1e+06")
  }
  tol <- .ckn(tol, "tol", f)
  if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol < 1e-12 ||
        tol > 1) {
    .ckstop(f, "`tol` must be one number from 1e-12 to 1")
  }
  list(maxit = maxit, tol = tol)
}

# mcartst(): the numeric columns as a matrix, columns and rows with no
# value dropped
.mcmat <- function(data, f) {
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) {
    .ckstop(f, "`data` must be a data frame or a numeric matrix")
  }
  num <- vapply(data, is.numeric, logical(1))
  if (!all(num)) {
    warning(sprintf("%s: dropping non-numeric column(s): %s.", f,
                    toString(names(data)[!num])), call. = FALSE)
    data <- data[, num, drop = FALSE]
  }
  if (ncol(data) < 2L) {
    .ckstop(f, "Little's test needs at least two numeric columns")
  }
  X <- .nummat(data, names(data))
  colnames(X) <- names(data)
  X <- X[, !apply(is.na(X), 2L, all), drop = FALSE]
  X <- X[!apply(is.na(X), 1L, all), , drop = FALSE]
  if (any(is.infinite(X))) {
    .ckstop(
      f,
      paste0(
        "`data` has infinite values, which Little's test cannot use; set ",
        "them to NA first"
      )
    )
  }
  if (ncol(X) < 2L || nrow(X) < 2L) {
    .ckstop(f, "too few complete columns or rows remain for Little's test")
  }
  X
}

# mcartst(): the note on a test with nothing to test or skipped patterns
.mcnote <- function(o) {
  df <- o[2L]
  note <- NA_character_
  if (df <= 0) {
    note <- paste0(
      "df is ", df, ": with a single missingness pattern there is nothing ",
      "to test. A complete table trivially satisfies MCAR."
    )
  }
  if (o[8L] > 0) {
    note <- paste(stats::na.omit(c(note, sprintf(
      paste0(
        "%d pattern(s) were skipped: their observed variables were ",
        "collinear."
      ),
      as.integer(o[8L])
    ))), collapse = " ")
  }
  note
}

.profrow <- function(v, nm) {
  nmiss <- sum(is.na(v))
  base <- data.frame(column = nm, type = class(v)[1L], n = length(v),
                     n_missing = nmiss,
                     pct_missing = 100 * nmiss / max(1L, length(v)),
                     stringsAsFactors = FALSE)
  na <- list(n_distinct = 0L, nzero = NA_integer_, nneg = NA_integer_,
             ninf = NA_integer_, mean = NA_real_, sd = NA_real_,
             median = NA_real_, mad = NA_real_, min = NA_real_,
             max = NA_real_, skewness = NA_real_, n_outliers = NA_integer_,
             hist = NA_character_, top = NA_character_)
  if (!is.numeric(v)) {
    ok <- as.character(v[!is.na(v)])
    na$n_distinct <- length(unique(ok))
    if (length(ok)) na$top <- names(sort(table(ok), decreasing = TRUE))[1L]
    return(cbind(base, as.data.frame(na, stringsAsFactors = FALSE)))
  }
  ok <- v[!is.na(v)]
  if (!length(ok)) {
    return(cbind(base, as.data.frame(na, stringsAsFactors = FALSE)))
  }
  fin <- as.numeric(ok[is.finite(ok)])
  na$n_distinct <- length(unique(ok))
  na$nzero <- sum(ok == 0)
  na$nneg <- sum(ok < 0)
  na$ninf <- sum(is.infinite(v))
  if (length(fin)) {
    s <- .Call(C_clab_r_pstat, fin)
    na$mean <- s[1L]
    na$sd <- s[2L]
    na$skewness <- s[3L]
    na$median <- s[4L]
    na$mad <- s[5L]
    na$min <- min(fin)
    na$max <- max(fin)
    na$n_outliers <- sum(fin < s[6L] | fin > s[7L])
    na$hist <- .spark(fin, 10L)
  }
  cbind(base, as.data.frame(na, stringsAsFactors = FALSE))
}

# dstntst(): a column of `data`, named by one string
.dtcol <- function(x, what, data, f) {
  .cktxt(x, what, f)
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !nzchar(trimws(x))) {
    .ckstop(f, "`%s` must be a single non-empty string", what)
  }
  if (!x %in% names(data)) .ckstop(f, "column `%s` not found in `data`", x)
  x
}

# dstntst(): `pooled`, NULL or labels of groups of `data`
.dtpool <- function(pooled, g, f) {
  if (is.null(pooled)) {
    return(invisible())
  }
  .cknofn(pooled, "pooled", f)
  .cktxt(pooled, "pooled", f)
  if (!is.character(pooled) || anyNA(pooled)) {
    .ckstop(
      f,
      "`pooled` must be the labels of the rows that pool other groups"
    )
  }
  miss <- setdiff(trimws(pooled), trimws(g))
  if (length(miss)) {
    .ckstop(f, "`pooled` names %s, not a group of `data`",
            paste0("\"", miss, "\"", collapse = ", "))
  }
}

# dstntst(): an estimate outside its own interval is not compared
.dtest <- function(y, lo, hi, g, f) {
  bad <- is.finite(lo) & is.finite(hi) & is.finite(y) & (lo > y | hi < y)
  if (any(bad)) {
    warning(sprintf(
      "%s: the interval does not contain the estimate for: %s; not compared.",
      f, toString(unique(g[bad]))
    ), call. = FALSE)
    y[bad] <- NA_real_
  }
  y
}

# dstntst(): the quality flags, blank as "ok" and an unknown one warned of
.dtqual <- function(q, f) {
  q <- trimws(as.character(q))
  q[is.na(q) | q == ""] <- "ok"
  odd <- setdiff(unique(q), c("ok", "E", "F", "x", "..", "..."))
  if (length(odd)) {
    warning(sprintf(
      "%s: quality flag(s) %s not recognised: treated as ok.", f,
      paste0("\"", odd, "\"", collapse = ", ")
    ), call. = FALSE)
    q[q %in% odd] <- "ok"
  }
  q
}

# dstntst(): the rows that pool other groups (`pooled`, or by their label)
.dtpan <- function(g, pooled) {
  if (!is.null(pooled)) {
    trimws(g) %in% trimws(pooled)
  } else {
    grepl(
      paste0(
        "^total|^indigenous identity|^aboriginal identity|^total ",
        "indigenous|^indigenous population"
      ),
      trimws(tolower(g))
    )
  }
}
