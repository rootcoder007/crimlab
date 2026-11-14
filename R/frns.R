#' Genotype match probability
#'
#' The chance that an unrelated person from the same population shares a DNA
#' profile, locus by locus and overall. With `theta = 0` this is the product
#' rule (`p^2` for a homozygote, `2pq` for a heterozygote). A positive `theta`
#' applies the Balding-Nichols correction for shared ancestry. That never
#' makes a homozygote match rarer than the product rule says; it can make a
#' heterozygote match rarer when both alleles are common. See
#' `proof("gtmp")`.
#'
#' @param p Frequency of the first allele at each locus: numbers above 0 and
#'   at most 1, optionally named by locus, as a numeric vector, a 1-d table, a
#'   list of single numbers or a one-row data frame.
#' @param q Frequency of the second allele at each locus; not read where
#'   `hom` is `TRUE`. Named inputs are matched to `p` by locus name, unnamed
#'   ones by position (only names protect against a mismatch).
#' @param hom `TRUE` where the genotype is homozygous: one value, or one per
#'   locus, matched to `p` by name when both are named.
#' @param theta The coancestry coefficient: one number from 0 up to, not
#'   including, 1.
#' @return A list: `locus`, a data frame with `locus`, `hom` and `mp`;
#'   `log10mp`, the base-10 logarithm of the profile match probability; and
#'   `log10lr`, the same with the sign changed, the weight of a single-source
#'   match against an unrelated alternative.
#' @references Balding, D. J. and Nichols, R. A. (1994). DNA profile match
#'   probability calculation: how to allow for population stratification,
#'   relatedness, database selection and single bands. *Forensic Science
#'   International*, 64(2-3),
#'   125-140. \doi{10.1016/0379-0738(94)90222-4}
#' @examples
#' gtmp(p = c(D3 = 0.1, vWA = 0.2), q = c(0.3, NA), hom = c(FALSE, TRUE),
#'      theta = 0.03)
#' @export
gtmp <- function(p, q = NULL, hom = FALSE, theta = 0) {
  .Call(C_clab_r_gtmp, p, q, hom, theta)
}

#' Log-likelihood-ratio cost
#'
#' Scores a set of likelihood ratios against the truth. Ratios from
#' same-source comparisons should be large and ratios from different-source
#' comparisons small. The cost is 1 for a system that always answers 1. It
#' is above 0 for any set of ratios, and tends to 0 as the ratios of a
#' system that separates the two kinds perfectly grow without bound. See
#' `proof("cllr")` and [cllrmin()].
#'
#' @param same Likelihood ratios from same-source comparisons: positive
#'   numbers.
#' @param diff Likelihood ratios from different-source comparisons: positive
#'   numbers.
#' @return One number, the cost in bits.
#' @references Brummer, N. and du Preez, J. (2006). Application-independent
#'   evaluation of speaker detection. *Computer Speech and Language*, 20(2-3),
#'   230-275. \doi{10.1016/j.csl.2005.08.001}
#' @examples
#' cllr(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#' @export
cllr <- function(same, diff) {
  .Call(C_clab_r_cllr, same, diff)
}

#' Posterior probability from a likelihood ratio
#'
#' Turns a likelihood ratio into the probability of the first hypothesis,
#' which needs a prior. A large ratio with a small prior can still leave that
#' probability low; reading the ratio as the probability is the prosecutor's
#' fallacy. See `proof("lrpost")`.
#'
#' @param lr Likelihood ratios: positive numbers.
#' @param prior The prior probability of the first hypothesis, strictly
#'   between 0 and 1: one value, or one per ratio, matched to `lr` by name
#'   when both are named with unique names.
#' @return The posterior probabilities.
#' @examples
#' lrpost(lr = 1e6, prior = 1 / 5e6)
#' @export
lrpost <- function(lr, prior) {
  .Call(C_clab_r_lrpost, lr, prior)
}
