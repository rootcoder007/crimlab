#' Summary statistics accumulated in blocks
#'
#' An accumulator for the count, mean, variance, skewness and kurtosis of a
#' numeric column that does not fit in memory, or that arrives in pieces (a
#' file read in chunks, results from several machines). Start one with
#' \code{onlsmry()}, feed it blocks with \code{smryupd()}, combine two with
#' \code{smrymrg()} and read the result with \code{smrystat()}.
#'
#' @details
#' The merge is exact algebra, not an approximation: it uses Chan, Golub and
#' LeVeque's parallel combination of central sums, extended to the third
#' and fourth moments by Terriberry, so accumulating a column in blocks
#' gives the answer of one pass over the whole column, to the rounding of
#' double arithmetic. Each block's sums are taken in two passes, as
#' \code{var()} takes them, so values far from zero keep their digits.
#' Memory stays constant in the number of blocks.
#'
#' The state is \eqn{(n, \bar{x}, M_2, M_3, M_4)} with
#' \eqn{M_k = \sum (x_i - \bar{x})^k}. \code{smrystat()} reports the sample
#' variance \eqn{M_2 / (n - 1)}, the moment skewness
#' \eqn{g_1 = (M_3/n) / (M_2/n)^{3/2}} and the excess kurtosis
#' \eqn{g_2 = (M_4/n) / (M_2/n)^2 - 3} (the biased moment estimators,
#' without small-sample correction). \code{NA} values are dropped from every
#' block. Undefined quantities are \code{NaN}: the mean for \eqn{n = 0}, the
#' variance for \eqn{n < 2}, the skewness for \eqn{n < 3} and the kurtosis
#' for \eqn{n < 4}, and both shape statistics for a constant column.
#'
#' The object is a plain list and is copied on assignment like any other R
#' value, so \code{smryupd()} returns the new accumulator, which must be
#' kept: \code{acc <- smryupd(acc, block)}.
#'
#' Errors: a block that is not numeric or holds an infinite value; an
#' accumulator that does not come from \code{onlsmry()}.
#'
#' @param x Numeric vector (finite; \code{NA} is skipped). Text that reads
#'   as numbers is taken as those numbers.
#' @param acc,a,b Accumulators from \code{onlsmry()}.
#' @inheritParams freqtbl
#' @return \code{onlsmry()}, \code{smryupd()} and \code{smrymrg()}: an
#'   accumulator, a list of class \code{onlsmry} with \code{state}.
#'   \code{smrystat()}: a named numeric vector \code{n}, \code{mean},
#'   \code{variance}, \code{sd}, \code{skewness}, \code{kurtosis}.
#' @references Chan, T. F., Golub, G. H. and LeVeque, R. J. (1983).
#'   Algorithms for computing the sample variance: analysis and
#'   recommendations. The American Statistician 37(3), 242-247.
#'   \doi{10.1080/00031305.1983.10483115}
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(1000) # simulated
#' # accumulate in blocks of 100
#' acc <- onlsmry()
#' for (i in seq(1, 1000, by = 100)) acc <- smryupd(acc, x[i:(i + 99)])
#' smrystat(acc)
#' all.equal(smrystat(acc)[["variance"]], stats::var(x))
#' # two accumulators built independently can be merged
#' a <- smryupd(onlsmry(), x[1:400])
#' b <- smryupd(onlsmry(), x[401:1000])
#' all.equal(smrystat(smrymrg(a, b)), smrystat(acc))
#' # an empty accumulator reports nothing rather than zero
#' smrystat(onlsmry())
#' try(onlsmry(c("a", "b")))
#' @export
onlsmry <- function(x = NULL) {
  acc <- structure(list(state = c(0, 0, 0, 0, 0)), class = c("onlsmry", "list"))
  if (is.null(x)) return(acc)
  .skupd(acc, .sknum(x, "x", "onlsmry"), "onlsmry")
}

#' @rdname onlsmry
#' @export
smryupd <- function(acc, x) {
  if (!inherits(acc, "onlsmry")) {
    .ckstop("smryupd", "`acc` must come from onlsmry()")
  }
  .skupd(acc, .sknum(x, "x", "smryupd"), "smryupd")
}

#' @rdname onlsmry
#' @export
smrymrg <- function(a, b) {
  if (!inherits(a, "onlsmry") || !inherits(b, "onlsmry")) {
    .ckstop("smrymrg", "both arguments must come from onlsmry()")
  }
  a$state <- .Call(C_clab_r_mommrg, a$state, b$state)
  a
}

#' @rdname onlsmry
#' @export
smrystat <- function(acc) {
  if (!inherits(acc, "onlsmry")) {
    .ckstop("smrystat", "`acc` must come from onlsmry()")
  }
  stats::setNames(
    .Call(C_clab_r_momstat, acc$state),
    c("n", "mean", "variance", "sd", "skewness", "kurtosis")
  )
}

#' @rdname onlsmry
#' @export
print.onlsmry <- function(x, ...) {
  s <- smrystat(x)
  cat(
    "Accumulated summary: n ", format(s[["n"]], big.mark = ","),
    ", mean ", formatC(s[["mean"]], format = "g", digits = 6),
    ", sd ", formatC(s[["sd"]], format = "g", digits = 6),
    ", variance ", formatC(s[["variance"]], format = "g", digits = 6),
    ", skewness ", formatC(s[["skewness"]], format = "g", digits = 4),
    ", excess kurtosis ", formatC(s[["kurtosis"]], format = "g", digits = 4),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' Uniform sample of a stream in one pass
#'
#' Vitter's Algorithm R: draws \code{k} items from a stream, each item
#' equally likely to be retained, in a single pass and using memory
#' proportional to \code{k} rather than to the stream's length.
#' \code{rsvidx()} returns the retained positions, so the same sample can be
#' applied to a file read line by line; \code{rsvsmp()} applies it to a
#' vector in memory.
#'
#' @details
#' The first \code{k} items fill the reservoir; item \eqn{i > k} then
#' replaces a uniformly chosen slot with probability \eqn{k / i}, so after
#' \eqn{n} items every item is in the sample with probability \eqn{k / n}
#' (Vitter 1985). The positions are returned sorted.
#'
#' The stream is sampled with the package's own generator seeded by
#' \code{seed}, not R's, so the sample is reproducible and R's random
#' stream is left alone. The generator is C++ \code{std::mt19937_64}, whose
#' output sequence the C++ standard fixes, and each draw is Lemire's
#' unbiased multiply-and-reject on its 64-bit outputs (Lemire 2019), so the
#' same seed gives the same sample on every platform and C++ library.
#'
#' Errors: \code{n} or \code{k} not one non-negative whole number; \code{n}
#' or \code{k} above 2147483647; \code{seed} not a whole number from 0 to
#' below \eqn{2^{53}}.
#'
#' @param x Vector to sample from.
#' @param n Length of the stream, one non-negative number.
#' @param k Sample size, one non-negative number. Values above the stream's
#'   length return the whole stream.
#' @param seed Seed for the generator, a whole number from 0 to below
#'   \eqn{2^{53}}.
#' @return \code{rsvidx()} a sorted numeric vector of positions (from 1);
#'   \code{rsvsmp()} the elements of \code{x} at them.
#' @references Vitter, J. S. (1985). Random sampling with a reservoir. ACM
#'   Transactions on Mathematical Software 11(1), 37-57.
#'   \doi{10.1145/3147.3165}
#'
#'   Lemire, D. (2019). Fast random integer generation in an interval. ACM
#'   Transactions on Modeling and Computer Simulation 29(1), 3.
#'   \doi{10.1145/3230636}
#' @examples
#' rsvidx(1000, 5, seed = 1)
#' rsvsmp(letters, 4, seed = 2)
#' rsvidx(3, 10)
#' # sample rows of a data frame
#' df <- data.frame(id = 1:20, v = (1:20)^2)
#' df[rsvidx(nrow(df), 3, seed = 7), ]
#' try(rsvidx(10, -1))
#' @export
rsvidx <- function(n, k, seed = 42L) {
  f <- "rsvidx"
  n <- .sknum(n, "n", f)
  k <- .sknum(k, "k", f)
  .rswhole(n, "n", f)
  .rswhole(k, "k", f)
  seed <- .rsseed(.sknum(seed, "seed", f), f)
  .Call(C_clab_r_rsv, n, k, seed)
}

#' @rdname rsvidx
#' @export
rsvsmp <- function(x, k, seed = 42L) {
  x[rsvidx(length(x), k, seed)]
}

#' Distinct-value count in fixed memory
#'
#' HyperLogLog: estimates how many distinct values a column holds using a
#' fixed register per bucket, regardless of the column's length or
#' cardinality. \code{dstnskch()} builds or extends the registers,
#' \code{dstncnt()} turns them into an estimate, and \code{skchmrg()}
#' combines two sketches.
#'
#' @details
#' Each value is hashed and the first \code{p} bits pick one of
#' \eqn{m = 2^p} registers; the register keeps the longest run of leading
#' zeros (plus one) seen in the remaining bits. The estimate is the
#' bias-corrected harmonic mean \eqn{\alpha_m m^2 / \sum_j 2^{-M_j}} of
#' Flajolet et al. (2007).
#'
#' The estimate carries a relative standard error of about
#' \code{1.04 / sqrt(2^p)}, so 0.8% at the default. It is an estimate:
#' \code{length(unique(x))} is exact when the column fits in memory. Below
#' roughly \code{2.5 * 2^p} distinct values the estimator switches to linear
#' counting, which is near-exact in that range.
#'
#' Values are hashed with SHA-256 (FIPS 180-4) of their UTF-8 text, so the
#' sketch is identical on every platform and across sessions. Values are
#' compared as text (\code{as.character()}), so \code{1}, \code{1.0} and
#' \code{"1"} are one value; \code{NA} is skipped. Since a register only
#' keeps a maximum, folding chunks with \code{regstr} and merging separate
#' sketches with \code{skchmrg()} give the same registers as one sketch of
#' all the data.
#'
#' Errors: \code{p} not a whole number from 4 to 20; \code{regstr} not of
#' length \eqn{2^p}; registers that are not whole numbers from 0 to 64; two
#' sketches of different sizes.
#'
#' @param x A vector, compared as text.
#' @param p Log2 of the register count, a whole number from 4 to 20
#'   (default 14).
#' @param regstr Registers from a previous call (with the same \code{p}), to
#'   extend; \code{NULL} starts afresh. For \code{dstncnt()}, the registers.
#' @param a,b Register sets to merge.
#' @return \code{dstnskch()} and \code{skchmrg()} an integer vector of
#'   \eqn{2^p} registers; \code{dstncnt()} the estimated number of distinct
#'   values.
#' @references Flajolet, P., Fusy, E., Gandouet, O. and Meunier, F. (2007).
#'   HyperLogLog: the analysis of a near-optimal cardinality estimation
#'   algorithm. Analysis of Algorithms 2007, 137-156.
#' @examples
#' set.seed(1)
#' x <- sample(1:5000, 200000, replace = TRUE) # simulated
#' dstncnt(dstnskch(x))
#' length(unique(x))
#' dstncnt(dstnskch(c("a", "b", "c", "a", "b")))
#' a <- dstnskch(x[1:100000])
#' b <- dstnskch(x[100001:200000])
#' dstncnt(skchmrg(a, b))
#' try(skchmrg(dstnskch(1:10, p = 4), dstnskch(1:10, p = 5)))
#' @export
dstnskch <- function(x, p = 14L, regstr = NULL) {
  f <- "dstnskch"
  p <- .dsp(.sknum(p, "p", f), f)
  if (!is.null(regstr)) regstr <- .dsreg(regstr, p, f)
  .Call(C_clab_r_hlladd, enc2utf8(as.character(x)), as.integer(p), regstr)
}

#' @rdname dstnskch
#' @export
dstncnt <- function(regstr) {
  if (is.null(regstr) || length(regstr) == 0L) return(0)
  .Call(C_clab_r_hllcnt, .skreg(regstr, "regstr", "dstncnt"))
}

#' @rdname dstnskch
#' @export
skchmrg <- function(a, b) {
  a <- .skreg(a, "a", "skchmrg")
  b <- .skreg(b, "b", "skchmrg")
  if (length(a) != length(b)) {
    .ckstop(
      "skchmrg", "both sketches must have the same number of registers"
    )
  }
  pmax(a, b)
}

# .rmbl_num(): numbers, or text that reads as numbers; not flags.
.sknum <- function(x, what, fn) {
  if (is.factor(x) || is.list(x) ||
        !(is.numeric(x) || is.logical(x) || is.character(x))) {
    .ckstop(fn, "`%s` must be numeric, got %s", what, class(x)[1L])
  }
  if (is.character(x)) {
    v <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & is.na(v) & nzchar(trimws(x)) & x != "NA"
    if (any(bad)) {
      .ckstop(
        fn, "`%s` must be numeric, not text such as \"%s\"", what, x[bad][1L]
      )
    }
    x <- v
  }
  if (is.logical(x) && !all(is.na(x))) {
    .ckstop(fn, "`%s` must be a number, not %s", what, format(x[!is.na(x)][1L]))
  }
  as.numeric(x)
}

.skupd <- function(acc, x, fn) {
  x <- x[!is.na(x)]
  if (any(is.infinite(x))) {
    .ckstop(
      fn,
      "`x` has infinite values, which no moment can hold; set them to NA first"
    )
  }
  if (!length(x)) return(acc)
  acc$state <- .Call(C_clab_r_momupd, acc$state, x)
  acc
}

# rsvidx(): `n` and `k`, whole numbers from 0 to 2147483647
.rswhole <- function(v, what, f) {
  if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 0 ||
        v != floor(v)) {
    .ckstop(f, "`%s` must be a single non-negative whole number", what)
  }
  if (v > .Machine$integer.max) {
    .ckstop(f, "`%s` must be at most 2147483647", what)
  }
}

.rsseed <- function(seed, f) {
  if (length(seed) != 1L || is.na(seed) || seed < 0 || seed != floor(seed) ||
        seed >= 2^53) {
    .ckstop(f, "`seed` must be a whole number in [0, 2^53)")
  }
  seed
}

.dsp <- function(p, f) {
  if (length(p) != 1L || is.na(p) || p < 4 || p > 20 || p != round(p)) {
    .ckstop(f, "`p` must be a single integer between 4 and 20")
  }
  p
}

# dstnskch(): the registers of a sketch with the same `p`
.dsreg <- function(regstr, p, f) {
  regstr <- .skreg(regstr, "regstr", f)
  if (length(regstr) != 2^p) {
    .ckstop(
      f,
      "`regstr` must have 2^p entries, from a sketch with the same `p`"
    )
  }
  regstr
}

.skreg <- function(r, what, fn) {
  if (is.factor(r) || is.list(r) || !(is.numeric(r) || is.logical(r))) {
    .ckstop(fn, "`%s` must be numeric, got %s", what, class(r)[1L])
  }
  r <- as.numeric(r)
  if (!all(is.finite(r)) || any(r < 0) || any(r > 64) || any(r != floor(r))) {
    .ckstop(fn, "`%s` must be non-negative integer registers", what)
  }
  as.integer(r)
}
