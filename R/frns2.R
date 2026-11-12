#' Kinship index
#'
#' The weight of two DNA profiles for a stated relationship against no
#' relationship, locus by locus and overall. Loci are taken as independent;
#' there is no allowance for mutation or coancestry. See `proof("kinidx")`.
#'
#' @param g1,g2 The two genotypes: two-column matrices of allele labels, one
#'   row per locus (row names name the loci), or two alleles for one locus.
#' @param freq The population frequencies of the alleles of `g1`, same shape.
#'   When `g1`, `g2` and `freq` all have row names, their rows are matched
#'   to `g1` by locus name (an error lists loci that differ); otherwise by
#'   position, so only names protect against a mismatch.
#' @param rel The relationship: `"parent"`, `"sibling"`, `"half"` (half
#'   siblings, grandparent or uncle) or `"cousin"`, or the three
#'   probabilities of sharing 0, 1 and 2 alleles by descent. Names are
#'   matched exactly, in lower case.
#' @return A list: `locus`, a data frame with `locus` and `idx`; `log10idx`,
#'   the base-10 logarithm of the combined index; and `post`, the probability
#'   of the relationship at even prior odds.
#' @examples
#' kinidx(g1 = c(12, 14), g2 = c(12, 14), freq = c(0.1, 0.2), rel = "sibling")
#' @export
kinidx <- function(g1, g2, freq, rel = "parent") {
  .Call(C_clab_r_kinidx, g1, g2, freq, rel)
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
#' @param mother `NULL`, or the mother's genotype in the same form. When all
#'   the tables have row names, their rows are matched to `child` by locus
#'   name; otherwise by position.
#' @return A list: `locus`, a data frame with `locus` and `idx` (`NaN` where
#'   the mother's genotype rules her out as mother); `log10idx`; and `post`,
#'   the probability of paternity at even prior odds.
#' @examples
#' patidx(child = c(12, 14), father = c(14, 15), freq = c(0.1, 0.2),
#'        mother = c(12, 12))
#' @export
patidx <- function(child, father, freq, mother = NULL) {
  .Call(C_clab_r_patidx, child, father, freq, mother)
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
  .Call(C_clab_r_mixlr, p, np, xp, nd, xd)
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
  .Call(C_clab_r_dbsrch, mp, n, pop)
}

#' Log-likelihood-ratio cost after recalibration
#'
#' Splits the cost of a set of likelihood ratios into `min`, their cost after
#' recalibration with the pool-adjacent-violators algorithm, which keeps
#' their order and gives each block of pooled ratios the value that costs
#' that block least, and `cal`, the rest (`cllr - min`). That `min` is the
#' least cost of every order-keeping recalibration, so that `cal` is never
#' negative, is not proved. A system can separate the two kinds of comparison
#' well and still report ratios that are too large or too small. See
#' `proof("cllrmin")` and `proof("lrcal")`.
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
  .Call(C_clab_r_cllrmin, same, diff)
}

#' Calibrated likelihood ratios
#'
#' Recalibrates a set of likelihood ratios against the truth with the
#' pool-adjacent-violators algorithm, keeping their order. The result can
#' hold 0 and `Inf`: the data never contradicted those ratios. See
#' `proof("lrcal")`.
#'
#' @inheritParams cllr
#' @return A list: `same` and `diff`, the calibrated ratios in the order
#'   given.
#' @examples
#' lrcal(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#' @export
lrcal <- function(same, diff) {
  .Call(C_clab_r_lrcal, same, diff)
}

#' Tippett plot
#'
#' For each value on the axis, the share of same-source and of
#' different-source likelihood ratios above it. A good system keeps the
#' same-source curve high to the right of 0 and the different-source curve
#' low; where each crosses 0 shows the rate of misleading evidence.
#' `proof()` records no theorem for this plot.
#'
#' @inheritParams cllr
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, a list of two data frames, `same` and `diff`, each with
#'   `log10lr` and `share` (the share of ratios at or above it).
#' @srrstats {EA5.0} Curves differ by line type, not colour, and no
#'   typeface or text size is set.
#' @srrstats {EA5.0b} No colour is used, so colour vision does not matter.
#' @examples
#' tipplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
#' @export
tipplt <- function(same, diff, ...) {
  r <- .Call(C_clab_r_tipplt, same, diff)
  s <- r[["same"]]
  d <- r[["diff"]]
  graphics::plot(range(s[["log10lr"]], d[["log10lr"]]), c(0, 1), type = "n",
                 xlab = "log10 likelihood ratio",
                 ylab = "share of ratios at or above", ...)
  graphics::lines(s[["log10lr"]], s[["share"]], type = "s", lty = 1)
  graphics::lines(d[["log10lr"]], d[["share"]], type = "s", lty = 2)
  graphics::abline(v = 0, lty = 3)
  graphics::legend("bottomleft", c("same source", "different source"),
                   lty = 1:2, bty = "n")
  invisible(r)
}

#' Empirical cross-entropy plot
#'
#' The information a set of likelihood ratios still leaves missing, across
#' prior odds: for the ratios as reported, for the same ratios after
#' calibration, and for a system that always reports 1. Reported ratios above
#' the last curve do harm at that prior. At even prior odds the first two
#' curves pass through the values of [cllrmin()]. `proof()` records no
#' theorem for this plot.
#'
#' @inheritParams cllr
#' @param lo The prior log10 odds to draw at.
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, a data frame: `lo`, `ece`, `min` and `ref`.
#' @references Ramos, D. and Gonzalez-Rodriguez, J. (2013). Reliable support:
#'   measuring calibration of likelihood ratios. *Forensic Science
#'   International*, 230(1-3), 156-169. \doi{10.1016/j.forsciint.2013.04.014}
#' @srrstats {EA5.0a} Text sizes are the device defaults; `cex` and other
#'   graphical parameters reach graphics::plot() through `...`.
#' @srrstats {EA5.4} Axis ticks come from graphics::plot(), via pretty().
#' @srrstats {EA5.5} The vertical axis is labelled in bits; the others are
#'   log10 odds and shares, which have no unit.
#' @examples
#' eceplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
#' @export
eceplt <- function(same, diff, lo = seq(-2.5, 2.5, by = 0.05), ...) {
  out <- .Call(C_clab_r_eceplt, same, diff, lo)
  graphics::plot(range(out[["lo"]]), c(0, max(out[["ece"]], out[["ref"]])),
                 type = "n", xlab = "prior log10 odds",
                 ylab = "empirical cross-entropy (bits)", ...)
  graphics::lines(out[["lo"]], out[["ece"]], lty = 1)
  graphics::lines(out[["lo"]], out[["min"]], lty = 2)
  graphics::lines(out[["lo"]], out[["ref"]], lty = 3)
  graphics::legend("topleft", c("reported", "calibrated", "always 1"),
                   lty = 1:3, bty = "n")
  invisible(out)
}
