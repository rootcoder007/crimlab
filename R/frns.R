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
#'   at most 1, optionally named by locus.
#' @param q Frequency of the second allele at each locus; not read where
#'   `hom` is `TRUE`.
#' @param hom `TRUE` where the genotype is homozygous: one value, or one per
#'   locus.
#' @param theta The coancestry coefficient: one number from 0 up to, not
#'   including, 1.
#' @return A list: `locus`, a data frame with `locus`, `hom` and `mp`;
#'   `log10mp`, the base-10 logarithm of the profile match probability; and
#'   `log10lr`, the same with the sign changed, the weight of a single-source
#'   match against an unrelated alternative.
#' @references Balding, D. J. and Nichols, R. A. (1994). DNA profile match
#'   probability calculation. *Forensic Science International*, 64(2-3),
#'   125-140. \doi{10.1016/0379-0738(94)90222-4}
#' @examples
#' gtmp(p = c(D3 = 0.1, vWA = 0.2), q = c(0.3, NA), hom = c(FALSE, TRUE),
#'      theta = 0.03)
#' @export
gtmp <- function(p, q = NULL, hom = FALSE, theta = 0) {
  nms <- .grp(p)
  p <- .num(p, "p", lo = 0, open = TRUE)
  n <- length(p)
  if (any(p > 1)) stop("`p` must be at most 1.", call. = FALSE)
  if (!is.logical(hom) || anyNA(hom) || !length(hom) %in% c(1L, n)) {
    stop("`hom` must be TRUE or FALSE, once or once per locus.", call. = FALSE)
  }
  hom <- rep_len(hom, n)
  if (is.null(q)) q <- rep(NA_real_, n)
  if (!is.numeric(q) && !all(is.na(q)) || length(q) != n) {
    stop("`q` must be numbers, one per locus.", call. = FALSE)
  }
  q <- as.double(q)
  q[hom] <- 0
  het <- !hom
  if (anyNA(q) || any(q[het] <= 0) || any(p[het] + q[het] > 1 + 1e-12)) {
    stop("`q` must be above 0 with `p + q` at most 1 at every heterozygous ",
         "locus.", call. = FALSE)
  }
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) ||
        theta < 0 || theta >= 1) {
    stop("`theta` must be one number from 0 up to, not including, 1.",
         call. = FALSE)
  }
  r <- .Call(C_csls_gtmp, p, q, hom, as.double(theta))
  list(locus = data.frame(locus = nms, hom = hom, mp = r[["locus"]]),
       log10mp = r[["log10mp"]], log10lr = -r[["log10mp"]])
}

#' Log-likelihood-ratio cost
#'
#' Scores a set of likelihood ratios against the truth. Ratios from
#' same-source comparisons should be large and ratios from different-source
#' comparisons small; the cost is 0 for a perfect system, 1 for one that
#' always answers 1, and above 1 for one that misleads. See `proof("cllr")`
#' and [cllrmin()].
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
  same <- .num(same, "same", lo = 0, open = TRUE)
  diff <- .num(diff, "diff", lo = 0, open = TRUE)
  .Call(C_csls_cllr, same, diff)
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
#'   between 0 and 1: one value, or one per ratio.
#' @return The posterior probabilities.
#' @examples
#' lrpost(lr = 1e6, prior = 1 / 5e6)
#' @export
lrpost <- function(lr, prior) {
  lr <- .num(lr, "lr", lo = 0, open = TRUE)
  prior <- .num(prior, "prior", lo = 0, open = TRUE)
  if (any(prior >= 1) || !length(prior) %in% c(1L, length(lr))) {
    stop("`prior` must be below 1, once or once per ratio.", call. = FALSE)
  }
  # 1 / (1 + 1 / odds): no overflow for a huge ratio
  1 / (1 + (1 - prior) / (prior * lr))
}
