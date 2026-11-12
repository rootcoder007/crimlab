#' Kinship index
#'
#' The weight of two DNA profiles for a stated relationship against no
#' relationship, locus by locus and overall. Loci are taken as independent;
#' there is no allowance for mutation or coancestry. See `proof("kinidx")`.
#'
#' @param g1,g2 The two genotypes: two-column matrices of allele labels, one
#'   row per locus (row names name the loci), or two alleles for one locus.
#' @param freq The population frequencies of the alleles of `g1`, same shape.
#' @param rel The relationship: `"parent"`, `"sibling"`, `"half"` (half
#'   siblings, grandparent or uncle) or `"cousin"`, or the three
#'   probabilities of sharing 0, 1 and 2 alleles by descent.
#' @return A list: `locus`, a data frame with `locus` and `idx`; `log10idx`,
#'   the base-10 logarithm of the combined index; and `post`, the probability
#'   of the relationship at even prior odds.
#' @examples
#' kinidx(g1 = c(12, 14), g2 = c(12, 14), freq = c(0.1, 0.2), rel = "sibling")
#' @export
kinidx <- function(g1, g2, freq, rel = "parent") {
  g1 <- .gt(g1, "g1")
  g2 <- .gt(g2, "g2", nrow(g1))
  f <- .gf(freq, g1, "g1")
  k <- .kin(rel)
  u <- unique(c(as.character(g1), as.character(g2)))
  .idx(g1, .Call(C_csls_kinidx, .code(g1, u), .code(g2, u), f, k))
}

#' Paternity index
#'
#' The weight of the DNA evidence that an alleged father is the father of a
#' child, against an unrelated man, locus by locus and overall. With the
#' mother's genotype the index uses the allele the child must have received
#' from the father; without it the index is the parent-child [kinidx()].
#' There is no allowance for mutation, so one inconsistent locus gives an
#' index of 0. See `proof("patidx")`.
#'
#' @param child,father Genotypes: two-column matrices of allele labels, one
#'   row per locus, or two alleles for one locus.
#' @param freq The population frequencies of the child's alleles, same shape.
#' @param mother `NULL`, or the mother's genotype in the same form.
#' @return A list: `locus`, a data frame with `locus` and `idx` (`NaN` where
#'   the mother's genotype rules her out as mother); `log10idx`; and `post`,
#'   the probability of paternity at even prior odds.
#' @examples
#' patidx(child = c(12, 14), father = c(14, 15), freq = c(0.1, 0.2),
#'        mother = c(12, 12))
#' @export
patidx <- function(child, father, freq, mother = NULL) {
  if (is.null(mother)) return(kinidx(child, father, freq, "parent"))
  child <- .gt(child, "child")
  n <- nrow(child)
  father <- .gt(father, "father", n)
  mother <- .gt(mother, "mother", n)
  f <- .gf(freq, child, "child")
  u <- unique(c(as.character(child), as.character(father),
                as.character(mother)))
  .idx(child, .Call(C_csls_patidx, .code(child, u), .code(mother, u),
                    .code(father, u), f))
}

#' Mixture likelihood ratio
#'
#' The weight of a DNA mixture at one locus for one account of who
#' contributed against another. Each account names the alleles of the
#' mixture its known contributors do not explain and the number of unknown
#' contributors who must explain them. Peak heights, drop-out and coancestry
#' are not modelled. See `proof("mixlr")`.
#'
#' With no alleles to explain and one unknown, the denominator is the chance
#' a random person is not excluded from the mixture.
#'
#' @param p The population frequencies of the alleles in the mixture: named
#'   numbers above 0 with a sum of at most 1.
#' @param np,nd The alleles left unexplained under the first and the second
#'   account: names of `p` (or `character(0)`), at most 20.
#' @param xp,xd The number of unknown contributors under each account.
#' @return A named numeric vector: `num`, `den` and `lr`.
#' @references Weir, B. S., Triggs, C. M., Starling, L., Stowell, L. I.,
#'   Walsh, K. A. J. and Buckleton, J. (1997). Interpreting DNA mixtures.
#'   *Journal of Forensic Sciences*, 42(2), 213-222.
#' @examples
#' # victim explains a and b; suspect would explain c and d
#' p <- c(a = 0.1, b = 0.2, c = 0.05, d = 0.15)
#' mixlr(p, np = character(0), xp = 0, nd = c("c", "d"), xd = 1)
#' @export
mixlr <- function(p, np, xp, nd, xd) {
  nms <- names(p)
  p <- .num(p, "p", lo = 0, open = TRUE)
  if (is.null(nms) || anyDuplicated(nms) || sum(p) > 1 + 1e-12) {
    stop("`p` must be uniquely named with a sum of at most 1.", call. = FALSE)
  }
  num <- .mix(p, nms, np, xp, "np", "xp")
  den <- .mix(p, nms, nd, xd, "nd", "xd")
  c(num = num, den = den, lr = num / den)
}

#' Database search
#'
#' What a match found by trawling a database is worth. `np` is the bound of
#' the second National Research Council report: the match probability times
#' the number of profiles searched. `any` is the exact chance of at least one
#' chance match among that many unrelated profiles, never above `np`. `post`
#' takes the other view (Balding and Donnelly): everyone else in the database
#' was excluded, so among `pop` possible sources with equal prior the matching
#' person's probability of being the source goes up, not down. See
#' `proof("dbsrch")`.
#'
#' @param mp The profile match probability: one number above 0 and at most 1.
#' @param n The number of profiles searched, the match included.
#' @param pop `NULL`, or the number of possible sources, at least `n`.
#' @return A named numeric vector: `np`, `any` and, when `pop` is given,
#'   `post`.
#' @references Balding, D. J. and Donnelly, P. (1996). Evaluating DNA profile
#'   evidence when the suspect is identified through a database search.
#'   *Journal of Forensic Sciences*, 41(4), 603-607.
#' @examples
#' dbsrch(mp = 1e-9, n = 5e6, pop = 3e7)
#' @export
dbsrch <- function(mp, n, pop = NULL) {
  mp <- .num(mp, "mp", len = 1L, lo = 0, open = TRUE)
  if (mp > 1) stop("`mp` must be at most 1.", call. = FALSE)
  n <- .num(n, "n", len = 1L, lo = 1)
  out <- c(np = min(1, n * mp), any = -expm1(n * log1p(-mp)))
  if (!is.null(pop)) {
    pop <- .num(pop, "pop", len = 1L, lo = n)
    out[["post"]] <- 1 / (1 + (pop - n) * mp)
  }
  out
}

#' Minimum log-likelihood-ratio cost
#'
#' Splits the cost of a set of likelihood ratios into what their ordering
#' alone would cost (`min`, found by recalibrating them with the
#' pool-adjacent-violators algorithm) and what their miscalibration adds
#' (`cal`). A system can separate the two kinds of comparison well and still
#' report ratios that are too large or too small.
#'
#' @inheritParams cllr
#' @return A named numeric vector: `cllr`, `min` and `cal`, with
#'   `cllr = min + cal`.
#' @references Brummer, N. and du Preez, J. (2006). Application-independent
#'   evaluation of speaker detection. *Computer Speech and Language*, 20(2-3),
#'   230-275. \doi{10.1016/j.csl.2005.08.001}
#' @examples
#' cllrmin(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#' @export
cllrmin <- function(same, diff) {
  r <- .cal(same, diff)
  full <- cllr(r[["same"]], r[["diff"]])
  c(cllr = full, min = r[["min"]], cal = full - r[["min"]])
}

#' Calibrated likelihood ratios
#'
#' Recalibrates a set of likelihood ratios against the truth with the
#' pool-adjacent-violators algorithm, keeping their order. The result can
#' hold 0 and `Inf`: the data never contradicted those ratios.
#'
#' @inheritParams cllr
#' @return A list: `same` and `diff`, the calibrated ratios in the order
#'   given.
#' @examples
#' lrcal(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#' @export
lrcal <- function(same, diff) {
  r <- .cal(same, diff)
  list(same = r[["csame"]], diff = r[["cdiff"]])
}

#' Tippett plot
#'
#' For each value on the axis, the share of same-source and of
#' different-source likelihood ratios above it. A good system keeps the
#' same-source curve high to the right of 0 and the different-source curve
#' low; where each crosses 0 shows the rate of misleading evidence.
#'
#' @inheritParams cllr
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, a list of two data frames, `same` and `diff`, each with
#'   `log10lr` and `share` (the share of ratios at or above it).
#' @examples
#' tipplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
#' @export
tipplt <- function(same, diff, ...) {
  same <- .num(same, "same", lo = 0, open = TRUE)
  diff <- .num(diff, "diff", lo = 0, open = TRUE)
  curve <- function(x) {
    x <- sort(log10(x))
    data.frame(log10lr = x, share = rev(seq_along(x)) / length(x))
  }
  s <- curve(same)
  d <- curve(diff)
  graphics::plot(range(s[["log10lr"]], d[["log10lr"]]), c(0, 1), type = "n",
                 xlab = "log10 likelihood ratio",
                 ylab = "share of ratios at or above", ...)
  graphics::lines(s[["log10lr"]], s[["share"]], type = "s", lty = 1)
  graphics::lines(d[["log10lr"]], d[["share"]], type = "s", lty = 2)
  graphics::abline(v = 0, lty = 3)
  graphics::legend("bottomleft", c("same source", "different source"),
                   lty = 1:2, bty = "n")
  invisible(list(same = s, diff = d))
}

#' Empirical cross-entropy plot
#'
#' The information a set of likelihood ratios still leaves missing, across
#' prior odds: for the ratios as reported, for the same ratios after
#' calibration, and for a system that always reports 1. Reported ratios above
#' the last curve do harm at that prior. At even prior odds the first two
#' curves pass through the values of [cllrmin()].
#'
#' @inheritParams cllr
#' @param lo The prior log10 odds to draw at.
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, a data frame: `lo`, `ece`, `min` and `ref`.
#' @references Ramos, D. and Gonzalez-Rodriguez, J. (2013). Reliable support:
#'   measuring calibration of likelihood ratios. *Forensic Science
#'   International*, 230(1-3), 156-169. \doi{10.1016/j.forsciint.2013.04.014}
#' @examples
#' eceplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
#' @export
eceplt <- function(same, diff, lo = seq(-2.5, 2.5, by = 0.05), ...) {
  r <- .cal(same, diff)
  lo <- .num(lo, "lo")
  ece <- function(s, d) {
    vapply(lo, function(l) {
      odds <- 10^l
      pr <- odds / (1 + odds)
      (pr * mean(.l1pinv(s * odds)) +
         (1 - pr) * mean(log1p(d * odds))) / log(2)
    }, numeric(1))
  }
  out <- data.frame(lo = lo, ece = ece(r[["same"]], r[["diff"]]),
                    min = ece(r[["csame"]], r[["cdiff"]]), ref = ece(1, 1))
  graphics::plot(range(lo), c(0, max(out[["ece"]], out[["ref"]])), type = "n",
                 xlab = "prior log10 odds", ylab = "empirical cross-entropy",
                 ...)
  graphics::lines(lo, out[["ece"]], lty = 1)
  graphics::lines(lo, out[["min"]], lty = 2)
  graphics::lines(lo, out[["ref"]], lty = 3)
  graphics::legend("topleft", c("reported", "calibrated", "always 1"),
                   lty = 1:3, bty = "n")
  invisible(out)
}

# A genotype table: one row per locus, two alleles. Returns the labels as
# characters, row by row.
.gt <- function(x, nm, n = NULL) {
  if (is.null(dim(x)) && length(x) == 2L) x <- matrix(x, 1L)
  if (!is.matrix(x) || ncol(x) != 2L || anyNA(x) ||
        !(is.numeric(x) || is.character(x))) {
    stop("`", nm, "` must be a two-column matrix of alleles, one row per ",
         "locus, with no missing values.", call. = FALSE)
  }
  if (!is.null(n) && nrow(x) != n) {
    stop("`", nm, "` must have ", n, " rows, not ", nrow(x), ".",
         call. = FALSE)
  }
  x
}

# Frequencies of the alleles of a genotype table, same shape.
.gf <- function(freq, g, nm) {
  if (is.null(dim(freq)) && length(freq) == 2L) freq <- matrix(freq, 1L)
  if (!is.matrix(freq) || !is.numeric(freq) || !identical(dim(freq), dim(g)) ||
        anyNA(freq) || any(freq <= 0) || any(freq > 1)) {
    stop("`freq` must hold the frequencies of the alleles of `", nm,
         "`, same shape, above 0 and at most 1.", call. = FALSE)
  }
  if (any(g[, 1L] == g[, 2L] & freq[, 1L] != freq[, 2L])) {
    stop("`freq` gives two frequencies for one allele of a homozygote.",
         call. = FALSE)
  }
  as.double(t(freq))
}

.code <- function(g, u) {
  as.integer(t(matrix(match(as.character(g), u), ncol = 2L)))
}

.idx <- function(g, idx) {
  l10 <- sum(log10(idx))
  list(locus = data.frame(locus = .grp(g[, 1L]), idx = idx), log10idx = l10,
       post = 1 / (1 + 10^(-l10)))
}

# The identity-by-descent probabilities of a relationship, by name or given.
.kin <- function(rel) {
  ks <- list(parent = c(0, 1, 0), sibling = c(0.25, 0.5, 0.25),
             half = c(0.5, 0.5, 0), cousin = c(0.75, 0.25, 0))
  if (is.character(rel) && length(rel) == 1L && rel %in% names(ks)) {
    return(ks[[rel]])
  }
  if (is.numeric(rel) && length(rel) == 3L && !anyNA(rel) &&
        all(rel >= 0) && abs(sum(rel) - 1) < 1e-9) {
    return(as.double(rel))
  }
  stop("`rel` must be one of ", toString(names(ks)),
       ", or three probabilities that sum to 1.", call. = FALSE)
}

# One side of a mixture ratio: the chance for alleles `need` and `x` unknowns.
.mix <- function(p, nms, need, x, a, b) {
  if (!is.character(need) || anyNA(need) || !all(need %in% nms) ||
        length(unique(need)) > 20L) {
    stop("`", a, "` must be names of `p`, at most 20.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 ||
        x != round(x) || x > 50) {
    stop("`", b, "` must be a whole number from 0 to 50.", call. = FALSE)
  }
  .Call(C_csls_mixp, p, nms %in% need, as.integer(x))
}

# Calibrated ratios, in the order given.
.cal <- function(same, diff) {
  same <- .num(same, "same", lo = 0, open = TRUE)
  diff <- .num(diff, "diff", lo = 0, open = TRUE)
  lr <- c(same, diff)
  is <- rep(c(TRUE, FALSE), c(length(same), length(diff)))
  o <- order(lr)
  r <- .Call(C_csls_lrcal, lr[o], is[o])
  post <- numeric(length(lr))
  post[o] <- r[["post"]]
  cal <- post / (1 - post) * length(diff) / length(same)
  list(same = same, diff = diff, csame = cal[is], cdiff = cal[!is],
       min = r[["min"]])
}


# log(1 + 1 / x) without forming 1 / x where it would overflow.
.l1pinv <- function(x) {
  small <- x < 1
  out <- log1p(1 / x)
  out[small] <- log1p(x[small]) - log(x[small])
  out
}
