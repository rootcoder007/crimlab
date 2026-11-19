#' Continuous DNA mixture likelihood ratio (gamma model)
#'
#' The weight of a DNA mixture for one account of who contributed against
#' another, from the peak heights, not only from which alleles were seen.
#' Each account names the reference profiles it holds to be contributors
#' (`hp`, the prosecution's, usually the person of interest and any
#' undisputed contributor; `hd`, the defence's); the rest of the `nc`
#' contributors are unknown people whose genotypes are summed over. The
#' model's parameters are fitted by maximum likelihood under each account
#' separately and the ratio of the two maxima reported, marker by marker
#' and over all markers.
#'
#' The model is the gamma model of peak heights (Cowell et al. 2015;
#' Bleka, Storvik and Gill 2016). Contributor `k` has mixture proportion
#' `m_k` (the proportions sum to 1) and carries 0, 1 or 2 copies `n_k(a)`
#' of allele `a`. The height of the peak at `a` is gamma with shape
#' `s(a) / omega^2` and scale `mu omega^2`, so its mean is `mu s(a)` and a
#' single heterozygous contributor's peaks have coefficient of variation
#' `omega`, where `s(a) = sum_k m_k n_k(a) d(a)`. With `deg = TRUE`,
#' `d(a) = beta^((size(a) - 125)` `/ 100)`, degradation by fragment size in
#' bases (`beta` fitted, above 0); otherwise `d(a) = 1`. With
#' `stut = TRUE`, each allele seen passes the share `xi` of its signal to
#' the position one repeat shorter (backward stutter): `s(a)` becomes
#' `(1 - xi) s(a) + xi s(a + 1)`, `xi` fitted. A peak lower than the
#' analytical threshold `at` is not seen, with the chance the gamma gives
#' (drop-out). A peak seen where `s(a) = 0` is a drop-in, of chance `pc`
#' times the allele's frequency times an exponential density of rate
#' `lambda` for its height above `at`; a marker without drop-in has chance
#' `1 - pc`. Unknown contributors' genotypes are drawn from the
#' frequencies with the coancestry correction `theta` (Balding and
#' Nichols), given the alleles of every profile in `ref`; markers are
#' independent. The alleles not seen at a marker, those of the reference
#' profiles among them, are taken together as one allele that does not
#' stutter and has the largest size listed at the marker; the positions
#' one repeat below the alleles seen are kept as positions no one carries.
#'
#' This is the model of EuroForMix (Bleka, Storvik and Gill 2016) with
#' backward stutter only and its frequencies not normalised (the rest of a
#' marker is the chance of an allele not listed), where the person of
#' interest is a known non-contributor under `hd` (every profile in `ref`
#' conditions the unknowns). Its tests compare the log likelihood with
#' EuroForMix's at given parameters, which agree to 1e-6 (EuroForMix keeps
#' its genotype weights in single precision), and the fitted maxima.
#'
#' Limits, plainly. The result is only as good as the model: the number of
#' contributors is the user's assumption, not estimated; forward and
#' double stutter, locus-specific amplification efficiencies, peak
#' height variance that differs by marker, and relatives among the
#' contributors are not modelled; the frequencies are used as given
#' (sampling error in them is ignored). Taking the alleles not seen as
#' one allele is exact without stutter and degradation; with them, an
#' allele not seen neither stutters nor is sized by its own fragment.
#' Stutter needs numeric labels (repeat numbers) for the alleles seen.
#' The fit is a numerical maximisation
#' (BFGS over transformed parameters from several starts) and can find a
#' local maximum; the sum over the unknowns grows as (number of
#' genotypes)^(unknowns) per marker and is refused above 2 million
#' combinations, so at most 5 contributors and in practice 3 unknowns.
#' This is not a validated casework system; its likelihoods are checked
#' against enumeration by hand and against EuroForMix, not against a
#' laboratory's validation data. Give `par` to evaluate both accounts at
#' given parameters instead
#' of fitting them.
#'
#' @param peaks The peaks seen: a data frame with columns `locus`,
#'   `allele` and `height`, every height at least `at`. A marker of `freq`
#'   with no peak counts too (everything dropped out).
#' @param freq The allele frequencies: a data frame with columns `locus`,
#'   `allele` and `freq` (above 0, summing to at most 1 at a marker; the
#'   rest is the chance of an allele not listed), and `size`, the fragment
#'   size in bases, needed when `deg = TRUE`. Every peak's allele must be
#'   listed.
#' @param ref The reference profiles: `NULL`, or a data frame with columns
#'   `id`, `locus`, `a1` and `a2`. Every profile conditions the unknowns'
#'   genotypes when `theta` is above 0.
#' @param hp,hd The ids of the profiles of `ref` that contribute under the
#'   first and the second account (`NULL`: none). Contributor `k` of
#'   an account is its `k`-th named profile, then its unknowns `U1`,
#'   `U2`, ...
#' @param nc The number of contributors: one whole number from 1 to 5.
#' @param at The analytical threshold: one number above 0.
#' @param pc The drop-in chance per marker: from 0, below 1.
#' @param lambda The rate of the exponential distribution of drop-in
#'   heights above `at`: above 0.
#' @param theta The coancestry coefficient: from 0, below 1.
#' @param stut,deg `TRUE` to model backward stutter and degradation.
#' @param par `NULL` to fit, or the parameters for both accounts: named
#'   numbers `mu`, `omega`, `beta` (with `deg`), `xi` (with `stut`) and
#'   `m1`, ..., `m<nc>` summing to 1.
#' @return An object of class `gmixlr`, a list: `lr` and `log10lr`; `fit`,
#'   a data frame with one row per account (`hyp`, `loglik` the natural
#'   log likelihood at its maximum, `mu`, `omega`, `beta`, `xi` (`NA` when
#'   not modelled), `conv` (`NA` with `par`) and `iter`); `mix`, the
#'   mixture proportions (`hyp`, `who`, `m`); and `locus`, the base-10 log
#'   ratio by marker at each account's parameters.
#' @references
#'   Cowell, R. G., Graversen, T., Lauritzen, S. L. and Mortera, J. (2015).
#'   Analysis of forensic DNA mixtures with artefacts. *Journal of the
#'   Royal Statistical Society Series C: Applied Statistics*, 64(1), 1-48.
#'   \doi{10.1111/rssc.12071}
#'
#'   Bleka, O., Storvik, G. and Gill, P. (2016). EuroForMix: An open source
#'   software based on a continuous model to evaluate STR DNA profiles from
#'   a mixture of contributors with artefacts. *Forensic Science
#'   International: Genetics*, 21, 35-44.
#'   \doi{10.1016/j.fsigen.2015.11.008}
#'
#'   Balding, D. J. and Nichols, R. A. (1994). DNA profile match
#'   probability calculation: how to allow for population stratification,
#'   relatedness, database selection and single bands. *Forensic Science
#'   International*, 64(2-3), 125-140. \doi{10.1016/0379-0738(94)90222-4}
#' @examples
#' # a two-person mixture at two markers: is S a contributor, with the
#' # victim V, or are V and an unknown person?
#' peaks <- data.frame(locus = rep(c("D3", "vWA"), c(4, 3)),
#'                     allele = c(14, 15, 16, 17, 16, 17, 18),
#'                     height = c(1210, 380, 1150, 420, 1400, 1960, 610))
#' freq <- data.frame(locus = rep(c("D3", "vWA"), c(5, 5)),
#'                    allele = c(13, 14, 15, 16, 17, 15, 16, 17, 18, 19),
#'                    freq = c(0.05, 0.1, 0.25, 0.25, 0.2,
#'                             0.1, 0.2, 0.3, 0.2, 0.1))
#' ref <- data.frame(id = rep(c("V", "S"), each = 2),
#'                   locus = rep(c("D3", "vWA"), 2),
#'                   a1 = c(14, 16, 15, 17), a2 = c(16, 17, 17, 18))
#' gmixlr(peaks, freq, ref, hp = c("V", "S"), nc = 2, hd = "V")
#' @export
gmixlr <- function(
  peaks,
  freq,
  ref = NULL,
  hp = NULL,
  nc,
  hd = NULL,
  at = 50,
  pc = 0.05,
  lambda = 0.01,
  theta = 0,
  stut = FALSE,
  deg = FALSE,
  par = NULL
) {
  .Call(
    C_clab_r_gmixlr,
    peaks,
    freq,
    ref,
    hp,
    nc,
    hd,
    at,
    pc,
    lambda,
    theta,
    stut,
    deg,
    par
  )
}

#' @export
print.gmixlr <- function(x, ...) {
  f <- x[["fit"]]
  cat(
    "gmixlr: ",
    nrow(x[["locus"]]),
    " markers; log10 LR ",
    format(x[["log10lr"]]),
    ", LR ",
    format(x[["lr"]]),
    "\n",
    sep = ""
  )
  print(f, row.names = FALSE)
  cat("mixture proportions:\n")
  print(x[["mix"]], row.names = FALSE)
  invisible(x)
}
