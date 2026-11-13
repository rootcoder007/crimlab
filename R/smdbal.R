#' Covariate balance and overlap
#'
#' The standardised mean difference and variance ratio of each covariate
#' between the treated and the controls, before and after weighting, with
#' a verdict by a threshold; and the overlap of the propensity between the
#' arms.
#'
#' The standardised mean difference is `(m1 - m0) / sqrt((s1^2 + s0^2) /
#' 2)` (Imbens and Rubin 2015): means weighted by `w` when given, the
#' standard deviations always those of the unweighted arms, as 'cobalt'
#' computes it by default, so that weighting moves only the numerator. The
#' variance ratio is `s1^2 / s0^2`, weighted when `w` is given (with the
#' unbiased correction for weights). `bal` is `TRUE` when every
#' difference (weighted, with `w`) is below `thr` in absolute value; 0.1 is
#' the usual rule.
#'
#' Overlap (Cole and Hernan 2008): the propensity's range in each arm, the
#' share of units outside the range both arms share (`out`), and the
#' largest inverse probability weight (`wmax`); `ovl` is `TRUE` when no unit
#' lies outside the common range. The propensity is the logistic
#' regression of `d` on `x` unless `ps` gives it.
#'
#' @inheritParams cbatt
#' @param w `NULL` or weights from 0, one per unit (from [cbatt()], say).
#' @param thr The threshold of the verdict: one positive number.
#' @return An object of class `smdbal`, a list: `tab`, one row per
#'   covariate (`var`, `m1`, `m0`, `smd`, `vr`, and with `w` also `m1w`,
#'   `m0w`, `smdw`, `vrw`); `overlap`, one row (`min1`, `max1`, `min0`,
#'   `max0`, `out`, `wmax`); `bal`; `ovl`; and `thr`.
#' @references
#'   Imbens, G. W. and Rubin, D. B. (2015). *Causal Inference for
#'   Statistics, Social, and Biomedical Sciences*. Cambridge University
#'   Press. \doi{10.1017/CBO9781139025751}
#'
#'   Cole, S. R. and Hernan, M. A. (2008). Constructing inverse probability
#'   weights for marginal structural models. *American Journal of
#'   Epidemiology*, 168(6), 656-664. \doi{10.1093/aje/kwn164}
#' @examples
#' set.seed(1)
#' x <- data.frame(age = rnorm(200, 40, 10), prior = rpois(200, 2))
#' d <- rbinom(200, 1, plogis((x$age - 40) / 10))
#' smdbal(x, d)
#' w <- cbatt(rnorm(200), d, x, est = "ebal", tgt = "att", reps = 0)$w
#' smdbal(x, d, w = w)$tab
#' @export
smdbal <- function(x, d, w = NULL, ps = NULL, thr = 0.1) {
  .Call(C_clab_r_smdbal, x, d, w, ps, thr)
}

#' @export
print.smdbal <- function(x, ...) {
  cat("smdbal: balanced at ", format(x[["thr"]]), ": ", x[["bal"]],
      "; overlap: ", x[["ovl"]], "\n", sep = "")
  print(x[["tab"]], row.names = FALSE)
  print(x[["overlap"]], row.names = FALSE)
  invisible(x)
}
