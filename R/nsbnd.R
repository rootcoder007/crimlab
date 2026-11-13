#' A true base rate behind a noisily recorded one
#'
#' A recorded outcome is often a proxy: a rearrest for a reoffence, a
#' recorded crime for one committed. With false-positive noise `a`, the
#' chance a case without the outcome is recorded with it, and
#' false-negative noise `b`, the chance a case with it is not, a true base
#' rate `p` is recorded as `pobs = p (1 - b) + (1 - p) a`. Given only the
#' boxes the noise rates lie
#' in, `nsbnd()` gives the interval of true rates consistent with each
#' recorded rate, and says for each pair of groups whether the order of
#' their true rates follows from the records whatever the noise.
#'
#' The true rate `(pobs - a) / (1 - a - b)` falls as `a` grows and rises as
#' `b` grows, so its ends are reached at corners of the boxes: the lower at
#' the largest `a` (capped at `pobs`) and smallest `b`, the upper at the
#' smallest `a` and largest `b` (capped at `1 - pobs`). Both ends are
#' attained, so no narrower interval follows from the boxes alone; with
#' known noise rates (each box one point) the interval is the one true
#' rate. A recorded rate that no admissible noise can produce (below the
#' smallest `a`, or above `1 - ` the smallest `b`) has no interval (`NaN`).
#'
#' Two groups' true rates are ordered whatever the noise exactly when their
#' intervals do not overlap; when they overlap there is admissible noise
#' under which either is the larger. `bd` is the largest upper end of the
#' `b` box, `a` held, at which the lower recorded rate's interval stays
#' below the other's: how much under-recording the order survives. The
#' boxes are the analyst's claim about how the records were made; the
#' bounds are only as good as that claim. See `proof("nsbnd")`.
#'
#' @param pobs Recorded rates: numbers from 0 to 1, named or not.
#' @param a The box of the false-positive noise rate: one number `a` (the
#'   box from 0 to `a`) or two, its ends; from 0.
#' @param b The box of the false-negative noise rate, in the same way. The
#'   upper ends of `a` and `b` must add to less than 1.
#' @return A list: `bnd`, a data frame with one row per rate: `rate` (its
#'   name, or its position), `pobs`, `lo`, `hi` and `width`; `pairs`, one
#'   row per pair `i` before `j`: `i`, `j`, `ord` (`"lower"`: `i`'s true
#'   rate is below `j`'s for every admissible noise; `"higher"`;
#'   `"undecided"`; `NA` when either has no interval) and `bd` (`NA` when
#'   the rates are equal or no box keeps them ordered).
#' @references Fogliato, R., Chouldechova, A. and G'Sell, M. (2020).
#'   Fairness evaluation in presence of biased noisy labels. *Proceedings
#'   of the 23rd International Conference on Artificial Intelligence and
#'   Statistics*, PMLR 108, 2325-2336. \doi{10.48550/arXiv.2003.13808}
#' @examples
#' # rearrest rates of two groups; up to 10% false positives, 20% missed
#' nsbnd(c(a = 0.35, b = 0.55), a = 0.10, b = 0.20)
#' # with up to 40% missed, the order is no longer identified
#' nsbnd(c(a = 0.35, b = 0.55), a = 0.05, b = 0.40)$pairs
#' # known noise rates give the true rate itself
#' nsbnd(0.4 * 0.8 + 0.6 * 0.1, a = c(0.1, 0.1), b = c(0.2, 0.2))$bnd
#' @export
nsbnd <- function(pobs, a, b) {
  .Call(C_clab_r_nsbnd, pobs, a, b)
}

#' The false positive rate the other rates imply
#'
#' Chouldechova's (2017) identity: in any confusion table, the false
#' positive rate is `p / (1 - p) (1 - ppv) / ppv (1 - fnr)`, from the base
#' rate `p`, the positive predictive value `ppv` and the false negative
#' rate `fnr`. So two groups whose scores are equally calibrated (the same
#' `ppv`) and miss equally (the same `fnr`) but whose base rates differ
#' must have different false positive rates: no score, however it is
#' built, can equalise all three. The identity is arithmetic, true of any
#' tool on any data; it says which equalities cannot hold together, not
#' which one should give way. See `proof("impfpr")`.
#'
#' @param p Base rates, strictly between 0 and 1.
#' @param ppv Positive predictive values, above 0 and at most 1.
#' @param fnr False negative rates, from 0 and below 1. The three have one
#'   length (or length 1); named ones are matched by name.
#' @return The false positive rates, one per triple.
#' @references Chouldechova, A. (2017). Fair prediction with disparate
#'   impact: a study of bias in recidivism prediction instruments. *Big
#'   Data*, 5(2), 153-163. \doi{10.1089/big.2016.0047}
#' @examples
#' # the same ppv and fnr, different base rates: different fpr
#' impfpr(p = c(a = 0.3, b = 0.5), ppv = 0.6, fnr = 0.25)
#' @export
impfpr <- function(p, ppv, fnr) {
  .Call(C_clab_r_impfpr, p, ppv, fnr)
}

#' Which pairs a noisy ranking orders
#'
#' People (or areas, or courts) ranked by estimated scores are ordered
#' correctly whatever the noise only when their estimates differ by more
#' than their two interval half-widths together: below that, admissible
#' errors can reverse them, above it none can. `rnkres()` says which pairs
#' that holds for, and the share it fails for; the resolution of the
#' ranking is twice the largest half-width. A league table whose intervals
#' are wide orders few of its pairs (Goldstein and Spiegelhalter 1996).
#' See `proof("rnkres")`.
#'
#' @param est The estimates: finite numbers, named or not, at least two.
#' @param hw The half-width of each estimate's interval: one number from 0,
#'   or one per estimate.
#' @return A list: `id`, a logical matrix, `TRUE` where the pair's order is
#'   identified (`NA` on the diagonal), with the estimates' names; `npair`,
#'   `nunid` (the pairs not identified), `share` (their share) and `res`.
#' @references Goldstein, H. and Spiegelhalter, D. J. (1996). League tables
#'   and their limitations: statistical issues in comparisons of
#'   institutional performance. *Journal of the Royal Statistical Society,
#'   Series A*, 159(3), 385-443. \doi{10.2307/2983325}
#' @examples
#' rnkres(c(a = 0.2, b = 0.35, c = 0.8), hw = c(0.1, 0.1, 0.05))
#' @export
rnkres <- function(est, hw) {
  .Call(C_clab_r_rnkres, est, hw)
}
