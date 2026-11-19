# The drift battery: has a table's distribution moved between two periods
# or two releases? A changed file hash says the bytes moved; these tests
# say whether the DISTRIBUTION did. Before computing a year-over-year
# change from a new release of an agency's published table, check the
# periods it shares with the last release against that release: a
# reclassified offence or a revised count shows up here first.
#
# Each result is a named number vector with the inputs kept in the
# attribute "lakner" (after the author of the method, as in
# rmoriebricklayer), printed with them, and plain numbers under `[`.

#' @export
print.drftm <- function(x, ...) {
  cat(class(x)[1L], ": ", sep = "")
  v <- unclass(x)
  attributes(v) <- list(names = names(x))
  cat(paste(names(v), format(v, digits = 4L), sep = " = ", collapse = ", "),
      "\n")
  i <- attr(x, "lakner")
  cat("inputs:", paste(names(i), unlist(i), sep = " = ", collapse = ", "), "\n")
  if (!is.null(attr(x, "method"))) cat("method:", attr(x, "method"), "\n")
  invisible(x)
}

#' @export
"[.drftm" <- function(x, ...) {
  v <- unclass(x)
  attributes(v) <- list(names = names(x))
  v[...]
}

#' Two-sample Kolmogorov-Smirnov test (C backend)
#'
#' Tests whether two numeric samples come from the same distribution, for
#' example a column of last year's release of a table and the same column
#' of this year's. It needs no assumption about the shape of either
#' distribution and reacts to any kind of change: a shift, a change in
#' spread or a change in shape.
#'
#' The statistic is the largest vertical gap between the two empirical
#' distribution functions, with the asymptotic two-sided p-value
#' \eqn{2\sum_k (-1)^{k-1} e^{-2k^2t^2}}, \eqn{t = \sqrt{n_{\mathrm{eff}}}D},
#' computed for \eqn{t < 1} by the theta-function form, which converges
#' there. It is accurate for moderate samples and conservative for small
#' ones; for an exact small-sample p-value use [stats::ks.test()]. Ties are
#' handled by comparing the two functions at each distinct value, so tied
#' data gives no warning.
#'
#' Missing values are dropped. Text, a factor or a list, an infinite
#' value, or no value left in either sample is refused.
#'
#' @param x,y Numeric vectors, the reference and the new sample; logical
#'   vectors are read as 0/1.
#' @return A named numeric vector of class `drftks`: `statistic` (the KS
#'   \eqn{D}), `p_value` and `neff` (\eqn{n_xn_y/(n_x+n_y)}), with the
#'   sample sizes in the attribute `lakner`.
#' @seealso [drftpsi()], [drftchisq()], [drfthomog()].
#' @examples
#' set.seed(1)
#' ref <- stats::rnorm(200)
#' drftks(ref, stats::rnorm(200))
#' drftks(ref, stats::rnorm(200, mean = 0.8))
#' drftks(ref, stats::rnorm(200, sd = 2.5))
#' drftks(ref, ref)[["statistic"]]
#' # simulated monthly counts of an offence in two releases
#' drftks(c(41, 38, 52, 47, 44, 50), c(41, 38, 52, 47, 61, 66))
#' @export
drftks <- function(x, y) {
  x <- .drfin(x, "x", "drftks")
  y <- .drfin(y, "y", "drftks")
  out <- .pcall(19L, x, y)
  names(out) <- c("statistic", "p_value", "neff")
  .drm(out, "drftks", list(nx = length(x), ny = length(y)))
}

#' Population stability index and Jensen-Shannon divergence (C backend)
#'
#' Two effect sizes of how far the distribution of a numeric column has
#' moved from a reference sample. Both samples are binned on the
#' reference's quantiles and the two sets of bin shares compared. Unlike
#' a p-value they do not grow with the sample size, so they say how large
#' a shift is rather than whether it is detectable.
#'
#' Binning: values are rounded to 12 significant digits (so a formatting
#' difference between two releases does not move tied values across an
#' edge), the `bins + 1` type-7 quantiles of `x` are computed, duplicates
#' removed, and the interior quantiles become right-closed cut points with
#' open ends. When the reference has fewer than three distinct quantiles,
#' every distinct quantile is a cut point instead.
#'
#' The population stability index is \eqn{\sum_i (p_i - q_i)\log(p_i/q_i)}
#' over the reference shares \eqn{p_i} and the new shares \eqn{q_i}: the
#' symmetrised Kullback-Leibler divergence. The conventional reading, from
#' credit-risk monitoring where it originates, is below 0.1 stable, 0.1 to
#' 0.25 worth a look, above 0.25 a material shift; these bands are
#' heuristics with no calibrated null distribution, and on a few hundred
#' rows binning noise alone can exceed 0.25. The Jensen-Shannon divergence
#' \eqn{\tfrac12 KL(p\|m) + \tfrac12 KL(q\|m)}, \eqn{m = (p+q)/2}, is
#' bounded by \eqn{\log 2} nats, so it compares across columns with
#' different numbers of bins, and is finite when a bin is empty on one
#' side. Empty bins are floored at `eps` for the index only.
#'
#' @inheritParams drftks
#' @param bins The number of quantile bins of the reference, one whole
#'   number from 2 to 1,000,000; fewer are used when quantiles tie.
#' @param eps The floor of an empty bin in the index, in (0, 1). It is part
#'   of the statistic: compare runs at one `eps`.
#' @return A named numeric vector of class `drftpsi`: `psi` and `jsdiv`,
#'   with the inputs in the attribute `lakner`.
#' @references Wu D, Olson DL (2010). Enterprise risk management: coping
#'   with model risk in a large bank. *Journal of the Operational Research
#'   Society* 61(2), 179--190. \doi{10.1057/jors.2008.144}
#'
#' Lin J (1991). Divergence measures based on the Shannon entropy. *IEEE
#'   Transactions on Information Theory* 37(1), 145--151.
#'   \doi{10.1109/18.61115}
#' @seealso [drftks()].
#' @examples
#' set.seed(2)
#' ref <- stats::rnorm(500)
#' drftpsi(ref, stats::rnorm(500))
#' drftpsi(ref, stats::rnorm(500, mean = 1))
#' drftpsi(ref, ref)
#' drftpsi(c(1, 1, 1), c(9, 9, 9))[["jsdiv"]] <= log(2)
#' drftpsi(ref, stats::rnorm(500, mean = 1), bins = 4, eps = 0.01)
#' try(drftpsi(ref, ref, bins = 1))
#' @export
drftpsi <- function(x, y, bins = 10L, eps = 1e-6) {
  x <- .drfin(x, "x", "drftpsi")
  y <- .drfin(y, "y", "drftpsi")
  .psargs(bins, eps)
  out <- .pcall(20L, x, y, c(bins, eps))
  names(out) <- c("psi", "jsdiv")
  .drm(out, "drftpsi", list(nx = length(x), ny = length(y), bins = bins))
}

#' Chi-square test of a categorical column against a known distribution
#'
#' Tests whether the category counts of a new sample agree with a
#' reference distribution taken as known (proportions fixed by a
#' specification, or a published breakdown): Pearson's chi-square
#' goodness-of-fit test. When the reference is itself a sample, use
#' [drfthomog()], which carries the reference's own sampling error.
#'
#' Categories present in one argument and not the other are aligned by
#' name, so a vanished category registers as a shortfall. With \eqn{O_j}
#' the observed counts and \eqn{E_j = N e_j / \sum e} the reference scaled
#' to the observed total, \deqn{X^2 = \sum_j \frac{(O_j - E_j)^2}{E_j}}
#' over the categories with \eqn{E_j > 0}, on \eqn{k - 1} degrees of
#' freedom; the p-value is the upper tail, computed directly. A category
#' that occurs but which the reference gives probability zero contradicts
#' it outright: `statistic` is then `Inf` and `p_value` 0. With a single
#' category there is nothing to test: `df` is 0 and the p-values `NA`. The
#' approximation is poor below expected counts of about 5.
#'
#' @param obs Named counts of the new sample, or a factor or text to be
#'   tabulated.
#' @param expd Named reference counts or proportions, or a factor or text
#'   to be tabulated; rescaled to the total of `obs`.
#' @return A named numeric vector of class `drftchisq`: `statistic`, `df`,
#'   `p_value` and `logpvalue` (the log of the upper tail, which keeps its
#'   information where `p_value` underflows), with the inputs in the
#'   attribute `lakner`.
#' @seealso [drfthomog()], [stats::chisq.test()].
#' @examples
#' ref <- c(a = 50, b = 30, c = 20)
#' drftchisq(c(a = 100, b = 60, c = 40), ref)
#' drftchisq(c(a = 40, b = 60, c = 100), ref)
#' drftchisq(c("a", "a", "b", "b"), c("a", "b"))
#' drftchisq(c("a", "a", "b", "d"), c("a", "a", "b", "b"))
#' @export
drftchisq <- function(obs, expd) {
  o <- .drtab(obs, "obs", "drftchisq")
  e <- .drtab(expd, "expd", "drftchisq")
  lv <- union(names(o), names(e))
  o <- .dralign(o, lv)
  e <- .dralign(e, lv)
  if (sum(o) <= 0) stop("drftchisq: `obs` has no counts.", call. = FALSE)
  if (sum(e) <= 0) stop("drftchisq: `expd` has no counts.", call. = FALSE)
  r <- .pcall(21L, unname(o), unname(e))
  out <- if (r[3L] > 0) {
    c(statistic = Inf, df = max(1, r[2L]), p_value = 0, logpvalue = -Inf)
  } else if (r[2L] < 1) {
    c(statistic = r[1L], df = 0, p_value = NA_real_, logpvalue = NA_real_)
  } else {
    c(statistic = r[1L], df = r[2L],
      p_value = stats::pchisq(r[1L], r[2L], lower.tail = FALSE),
      logpvalue = stats::pchisq(r[1L], r[2L], lower.tail = FALSE,
                                log.p = TRUE))
  }
  .drm(out, "drftchisq", list(n = sum(o), k = length(lv)))
}

#' Chi-square test of homogeneity of two categorical samples
#'
#' Tests whether two samples of a categorical variable (offence type,
#' region, race in last year's and this year's release) were drawn from
#' the same distribution, by Pearson's chi-square on the 2-by-k table of
#' their counts.
#'
#' With pooled shares \eqn{\hat\pi_j = (x_j + y_j)/(n_x + n_y)},
#' \deqn{X^2 = \sum_j \frac{(x_j - n_x\hat\pi_j)^2}{n_x\hat\pi_j} +
#' \frac{(y_j - n_y\hat\pi_j)^2}{n_y\hat\pi_j}} on \eqn{k - 1} degrees of
#' freedom. Three branches, named in the attribute `method`:
#' `"chi-square approximation"` when every expected count is at least 5;
#' `"Monte Carlo permutation` `(small expected counts)"` otherwise, the
#' p-value from [stats::chisq.test()] with `simulate.p.value = TRUE` and
#' `B = 10000` tables with the same margins, drawn under `seed` with the
#' caller's random stream restored (`p_value` is then at least 1/(B+1));
#' and `"inapplicable: ..."` for a single category or fewer than two
#' observations per category, where no test has power and `p_value` is
#' `NA` rather than 1.
#'
#' [drftchisq()] takes the reference as known; when it is itself a finite
#' sample, that ignores its sampling error and reports drift too readily.
#' This test estimates the shared distribution from the pooled margins and
#' carries the uncertainty of both samples, the right test for two
#' releases of a table. Categories present in only one sample are given
#' zero counts in the other.
#'
#' @param x,y Factors or text (or named count vectors): the reference and
#'   the new sample.
#' @param seed The seed of the Monte Carlo branch, one whole number from
#'   -2147483647 to 2147483647.
#' @return A named numeric vector of class `drfthomog`: `statistic`, `df`,
#'   `p_value` and `logpvalue`, with the attribute `method` and the inputs
#'   in the attribute `lakner`.
#' @seealso [drftchisq()], [stats::chisq.test()].
#' @examples
#' set.seed(1)
#' a <- sample(c("x", "y", "z"), 300, TRUE)
#' b <- sample(c("x", "y", "z"), 300, TRUE)
#' drfthomog(a, b)
#' drfthomog(a, sample(c("x", "y", "z"), 300, TRUE, prob = c(0.7, 0.2, 0.1)))
#' drfthomog(c(a = 30, b = 20), c(a = 25, b = 25))
#' attr(drfthomog(c(a = 3, b = 1, c = 2), c(a = 1, b = 4, c = 1)), "method")
#' try(drfthomog(c(1, 2), c(a = 3)))
#' @export
drfthomog <- function(x, y, seed = 1L) {
  cx <- .drtab(x, "x", "drfthomog")
  cy <- .drtab(y, "y", "drfthomog")
  lv <- union(names(cx), names(cy))
  cx <- .dralign(cx, lv)
  cy <- .dralign(cy, lv)
  if (sum(cx) <= 0 || sum(cy) <= 0) {
    stop("drfthomog: both samples must have at least one observation.",
         call. = FALSE)
  }
  r <- .pcall(22L, unname(cx), unname(cy))
  st <- r[1L]
  df <- r[2L]
  n <- sum(cx) + sum(cy)
  inp <- list(nx = sum(cx), ny = sum(cy), k = length(lv))
  if (df < 1 || n / r[4L] < 2) {
    m <- if (df < 1) {
      "inapplicable: one category"
    } else {
      paste0("inapplicable: ", r[4L], " categories for ", n,
             " observations (fewer than two per category)")
    }
    out <- c(statistic = st, df = max(df, 0), p_value = NA_real_,
             logpvalue = NA_real_)
    return(.drm(out, "drfthomog", inp, m))
  }
  if (r[3L] < 5) {
    keep <- (cx + cy) > 0
    p <- .drseed(seed, suppressWarnings(stats::chisq.test(
      rbind(cx[keep], cy[keep]), simulate.p.value = TRUE, B = 10000L
    ))$p.value)
    out <- c(statistic = st, df = df, p_value = p, logpvalue = log(p))
    m <- "Monte Carlo permutation (small expected counts)"
  } else {
    out <- c(statistic = st, df = df,
             p_value = stats::pchisq(st, df, lower.tail = FALSE),
             logpvalue = stats::pchisq(st, df, lower.tail = FALSE,
                                       log.p = TRUE))
    m <- "chi-square approximation"
  }
  .drm(out, "drfthomog", inp, m)
}

#' A drift statistic `v` with its inputs and method as attributes, and its
#'   class.
#' @noRd
.drm <- function(v, name, inputs, method = NULL) {
  structure(v, lakner = inputs, method = method,
            class = c(name, "drftm", "numeric"))
}

#' `x` as finite numbers with the missing ones dropped; refused when infinite
#'   or empty.
#' @noRd
.drfin <- function(x, nm, fn) {
  x <- .pnum(x, nm, fn)
  x <- x[!is.na(x)]
  if (any(is.infinite(x))) {
    stop(fn, ": `", nm, "` must be finite (no Inf).", call. = FALSE)
  }
  if (!length(x)) {
    stop(fn, ": `", nm, "` needs at least one non-missing value, got 0.",
         call. = FALSE)
  }
  x
}

#' The arguments `bins` and `eps` of drftpsi() checked.
#' @noRd
.psargs <- function(bins, eps) {
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 2 ||
        bins != floor(bins)) {
    stop("drftpsi: `bins` must be one whole number of at least 2.",
         call. = FALSE)
  }
  if (bins > 1e6) {
    stop("drftpsi: `bins` must be at most 1,000,000.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0 ||
        eps >= 1) {
    stop("drftpsi: `eps` must be a positive floor below one (0 < eps < 1).",
         call. = FALSE)
  }
}

#' counts by label: a factor or text is tabulated, a numeric vector must be
#' named counts
#' @noRd
.drtab <- function(v, nm, fn) {
  if (is.factor(v) || is.character(v)) {
    t <- table(as.character(v))
    out <- stats::setNames(as.numeric(t), names(t))
  } else {
    if (is.null(names(v)) || !is.numeric(v)) {
      stop(fn, ": a numeric `", nm, "` must be named.", call. = FALSE)
    }
    out <- stats::setNames(as.numeric(v), names(v))
  }
  if (!all(is.finite(out)) || any(out < 0)) {
    stop(fn, ": `", nm, "` must be finite, non-negative counts.", call. = FALSE)
  }
  out
}

#' The counts `o` laid out over the levels `lv`, 0 where a level is absent.
#' @noRd
.dralign <- function(o, lv) stats::setNames(ifelse(is.na(o[lv]), 0, o[lv]), lv)

#' expr evaluated under set.seed(seed), the caller's stream restored
#' @noRd
.drseed <- function(seed, expr) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        seed != floor(seed) || abs(seed) > .Machine$integer.max) {
    stop("drfthomog: `seed` must be one whole number from -2147483647 to ",
         "2147483647.", call. = FALSE)
  }
  e <- globalenv()
  had <- exists(".Random.seed", envir = e, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = e, inherits = FALSE)
  on.exit(if (had) {
    assign(".Random.seed", old, envir = e)
  } else if (exists(".Random.seed", envir = e, inherits = FALSE)) {
    rm(".Random.seed", envir = e)
  }, add = TRUE)
  set.seed(as.integer(seed))
  expr
}
