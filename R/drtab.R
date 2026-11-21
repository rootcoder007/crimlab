#' Drift of a whole table between two releases
#'
#' Compares a reference table and a current one column by column, with the
#' tests of this family: a numeric (or date) column by the two-sample
#' Kolmogorov-Smirnov test ([drftks()]) and the population stability index
#' ([drftpsi()]), a text, factor or logical column by the chi-square test
#' of homogeneity ([drfthomog()]). `drftcal()` then says how often these
#' screens fire on a table that has not changed at all, which is what a
#' flag has to be read against.
#'
#' `drfttab()` screens the columns the two tables share. A column is
#' flagged when its test rejects at `alpha`, or when its population
#' stability index is above `psithr` and both samples have at least
#' `psiminn` values (on fewer, binning noise alone crosses the usual
#' bands). A column that cannot be tested (all missing on one side, one
#' category, of a type that is neither numbers nor labels) has `flag`
#' `NA`, and `note` says why. Missing values are left out of each test;
#' `nref` and `ncur` are the counts that entered it.
#'
#' Many columns mean many tests: at `alpha` 0.01, one column in a hundred
#' is flagged on unchanged data by the tests alone. `drftcal()` measures
#' this for the table at hand: it splits `data` at random into two halves
#' `n` times, runs `drfttab()` on each pair, and reports the share of
#' splits on which each column, and any column, was flagged. Each half has
#' half the rows, so the rates describe a comparison at that size.
#'
#' Both functions exist in every language of the package. One branch is
#' R's alone: a column of labels with an expected count below 5 gets its
#' p-value from a Monte Carlo draw of R's generator ([drfthomog()]); the
#' other languages leave such a column untested (`flag` missing) and say
#' so in `note`.
#'
#' Errors: `ref`, `cur` or `data` not a data frame; no shared column;
#' `alpha` outside (0, 1); `psithr` not a positive number; `psiminn`,
#' `bins` or `n` not a whole number of at least 1 (`bins` at least 2, `n`
#' at most 1,000,000); fewer than four rows in `data`; a `seed` that is
#' not one whole number (`drftcal()`: from 0).
#'
#' @param ref,cur The reference and the current table: data frames. Only
#'   the columns named in both are compared.
#' @param alpha The level of each column's test, in (0, 1).
#' @param psithr The population stability index above which a numeric
#'   column is flagged; 0.25 is the conventional "material shift".
#' @param psiminn The fewest values each sample needs for `psithr` to
#'   count.
#' @inheritParams drftpsi
#' @param seed The seed of the Monte Carlo p-value [drfthomog()] uses for
#'   small expected counts, and in `drftcal()` of the random splits too
#'   (there a whole number from 0: [cunif()] draws a seed for each split
#'   from it, and a split takes the rows of the smaller half of the
#'   uniforms drawn at its seed, so a seed gives the same splits in every
#'   language of the package).
#' @return `drfttab()` an object of class `drfttab`, a data frame with a
#'   row per shared column: `column`, `type` (`"numeric"` or `"label"`),
#'   `test`, `statistic`, `p`, `psi` (`NA` for a label column), `nref`,
#'   `ncur`, `flag` and `note`; the attributes `only` (the columns in one
#'   table alone), `alpha`, `psithr` and `psiminn`. `drftcal()` an object
#'   of class `drftcal`, a list: `columns`, a data frame of `column` and
#'   `far` (the share of splits on which it was flagged); `anyflag`, the
#'   share of splits with at least one flag; `alpha`; `n`; and `alphafam`,
#'   `alpha` divided by the number of columns, the level per column that
#'   holds the chance of any false flag at `alpha` (Bonferroni).
#' @references
#'   Dunn, O. J. (1961). Multiple comparisons among means. *Journal of the
#'   American Statistical Association*, 56(293), 52-64.
#'   \doi{10.1080/01621459.1961.10482090}
#' @seealso [drftks()], [drftpsi()], [drfthomog()] for one column.
#' @examples
#' set.seed(1)
#' ref <- data.frame(age = rnorm(300, 35, 9),
#'                   charge = sample(c("assault", "theft", "drugs"), 300, TRUE))
#' cur <- data.frame(age = rnorm(300, 39, 9),
#'                   charge = sample(c("assault", "theft", "drugs"), 300, TRUE))
#' drfttab(ref, cur)
#' # how often the screens fire when nothing changed
#' drftcal(ref, n = 20)
#' try(drfttab(ref, data.frame(other = 1)))
#' @export
drfttab <- function(
  ref,
  cur,
  alpha = 0.01,
  psithr = 0.25,
  psiminn = 100L,
  bins = 10L,
  seed = 1L
) {
  fn <- "drfttab"
  if (!is.data.frame(ref) || !is.data.frame(cur)) {
    .ckstop(fn, "`ref` and `cur` must be data frames")
  }
  a <- .drtargs(alpha, psithr, psiminn, fn)
  cols <- intersect(names(ref), names(cur))
  if (!length(cols)) {
    .ckstop(fn, "`ref` and `cur` share no column")
  }
  rows <- lapply(cols, function(k) {
    .clas(fn, .drtcol(ref[[k]], cur[[k]], a, bins, seed))
  })
  out <- data.frame(
    column = cols,
    do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)),
    stringsAsFactors = FALSE
  )
  structure(
    out,
    only = list(
      ref = setdiff(names(ref), cols),
      cur = setdiff(names(cur), cols)
    ),
    alpha = a[["alpha"]],
    psithr = a[["psithr"]],
    psiminn = a[["psiminn"]],
    class = c("drfttab", "data.frame")
  )
}

#' @rdname drfttab
#' @param data One table, known not to have drifted (one release).
#' @param n The number of random half-splits, from 1 to 1,000,000.
#' @export
drftcal <- function(
  data,
  n = 50L,
  alpha = 0.01,
  seed = 1L,
  psithr = 0.25,
  psiminn = 100L,
  bins = 10L
) {
  fn <- "drftcal"
  if (!is.data.frame(data) || nrow(data) < 4L) {
    .ckstop(fn, "`data` must be a data frame with at least four rows")
  }
  if (
    !is.numeric(n) ||
      length(n) != 1L ||
      !is.finite(n) ||
      n < 1 ||
      n > 1e6 ||
      n != floor(n)
  ) {
    .ckstop(fn, "`n` must be a whole number from 1 to 1000000")
  }
  if (
    !is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed < 0 ||
      seed != floor(seed) ||
      seed >= 2^53
  ) {
    .ckstop(fn, "`seed` must be one whole number from 0, below 2^53")
  }
  m <- nrow(data)
  # each split has a seed of its own, drawn from `seed` by the kernel's
  # generator; its half is the rows of the m %/% 2 smallest of m uniforms
  # at that seed: the same rows in every language
  each <- floor(cunif(n, seed) * 2^52)
  hit <- lapply(seq_len(n), function(i) {
    k <- sort(order(cunif(m, each[i]))[seq_len(m %/% 2L)])
    r <- .clas(
      fn,
      drfttab(
        data[k, , drop = FALSE],
        data[-k, , drop = FALSE],
        alpha = alpha,
        psithr = psithr,
        psiminn = psiminn,
        bins = bins,
        seed = seed
      )
    )
    stats::setNames(!is.na(r[["flag"]]) & r[["flag"]], r[["column"]])
  })
  h <- do.call(rbind, hit)
  structure(
    list(
      columns = data.frame(
        column = colnames(h),
        far = unname(colMeans(h)),
        stringsAsFactors = FALSE
      ),
      anyflag = mean(rowSums(h) > 0),
      alpha = alpha,
      n = as.integer(n),
      alphafam = alpha / ncol(h)
    ),
    class = "drftcal"
  )
}

#' @export
print.drfttab <- function(x, ...) {
  f <- x[["flag"]]
  cat(
    "drfttab: ",
    nrow(x),
    " columns, ",
    sum(f, na.rm = TRUE),
    " flagged, ",
    sum(is.na(f)),
    " not testable (alpha ",
    format(attr(x, "alpha")),
    ", PSI above ",
    format(attr(x, "psithr")),
    ")\n",
    sep = ""
  )
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  d[["statistic"]] <- signif(d[["statistic"]], 4L)
  d[["p"]] <- signif(d[["p"]], 3L)
  d[["psi"]] <- signif(d[["psi"]], 3L)
  print(d[, names(d) != "note"], row.names = FALSE)
  o <- attr(x, "only")
  if (length(o[["ref"]])) {
    cat("only in `ref`:", paste(o[["ref"]], collapse = ", "), "\n")
  }
  if (length(o[["cur"]])) {
    cat("only in `cur`:", paste(o[["cur"]], collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
print.drftcal <- function(x, ...) {
  cat(
    "drftcal: ",
    x[["n"]],
    " half-splits of unchanged data at alpha ",
    format(x[["alpha"]]),
    "\nany column flagged on ",
    format(round(100 * x[["anyflag"]], 1)),
    "% of splits; per-column level for a family-wise ",
    format(x[["alpha"]]),
    ": ",
    format(signif(x[["alphafam"]], 3L)),
    "\n",
    sep = ""
  )
  d <- x[["columns"]]
  d <- d[d[["far"]] > 0, , drop = FALSE]
  if (nrow(d)) {
    print(d[order(-d[["far"]]), , drop = FALSE], row.names = FALSE)
  } else {
    cat("no column was flagged on any split\n")
  }
  invisible(x)
}

#' alpha, psithr and psiminn of drfttab() checked, as a list.
#' @noRd
.drtargs <- function(alpha, psithr, psiminn, fn) {
  one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!one(alpha) || alpha <= 0 || alpha >= 1) {
    .ckstop(fn, "`alpha` must be one number in (0, 1)")
  }
  if (!one(psithr) || psithr <= 0) {
    .ckstop(fn, "`psithr` must be one positive number")
  }
  if (!one(psiminn) || psiminn < 1 || psiminn != floor(psiminn)) {
    .ckstop(fn, "`psiminn` must be a whole number of at least 1")
  }
  list(alpha = alpha, psithr = psithr, psiminn = psiminn)
}

#' One column of drfttab(): the row of the result for the reference values
#' x and the current values y.
#' @noRd
.drtcol <- function(x, y, a, bins, seed) {
  row <- function(type, test, st, p, psi, nx, ny, note = "") {
    flag <- if (is.na(p) && is.na(psi)) {
      NA
    } else {
      isTRUE(p < a[["alpha"]]) ||
        (isTRUE(psi > a[["psithr"]]) && min(nx, ny) >= a[["psiminn"]])
    }
    list(
      type = type,
      test = test,
      statistic = st,
      p = p,
      psi = psi,
      nref = nx,
      ncur = ny,
      flag = flag,
      note = note
    )
  }
  num <- function(v) {
    is.numeric(v) || inherits(v, c("Date", "POSIXt"))
  }
  lab <- function(v) is.character(v) || is.factor(v) || is.logical(v)
  # a column with nothing but missing values has no kind of its own: it
  # takes the other side's
  void <- function(v) is.atomic(v) && all(is.na(v))
  if (void(x) && !void(y)) {
    x <- y[rep(NA_integer_, length(x))]
  } else if (void(y) && !void(x)) {
    y <- x[rep(NA_integer_, length(y))]
  } else if (void(x) && void(y)) {
    x <- as.character(x)
    y <- as.character(y)
  }
  if (num(x) && num(y)) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (!length(x) || !length(y)) {
      return(row(
        "numeric",
        "ks",
        NA_real_,
        NA_real_,
        NA_real_,
        length(x),
        length(y),
        "no values on one side"
      ))
    }
    k <- drftks(x, y)
    return(row(
      "numeric",
      "ks",
      k[["statistic"]],
      k[["p_value"]],
      drftpsi(x, y, bins = bins)[["psi"]],
      length(x),
      length(y)
    ))
  }
  if (lab(x) && lab(y)) {
    x <- as.character(x)
    y <- as.character(y)
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if (!length(x) || !length(y)) {
      return(row(
        "label",
        "homogeneity",
        NA_real_,
        NA_real_,
        NA_real_,
        length(x),
        length(y),
        "no values on one side"
      ))
    }
    h <- drfthomog(x, y, seed = seed)
    return(row(
      "label",
      "homogeneity",
      h[["statistic"]],
      h[["p_value"]],
      NA_real_,
      length(x),
      length(y),
      if (is.na(h[["p_value"]])) attr(h, "method") else ""
    ))
  }
  row(
    "other",
    "none",
    NA_real_,
    NA_real_,
    NA_real_,
    length(x),
    length(y),
    "not numbers on both sides nor labels on both sides"
  )
}
