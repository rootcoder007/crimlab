#' Design-based causal effects of a 0/1 treatment
#'
#' The battery of estimators of the effect of a 0/1 treatment `d` on an
#' outcome `y` given covariates `x` that do not learn the nuisances by
#' cross-fitting: Hajek-normalised inverse probability weighting,
#' g-computation, subclassification on the propensity, nearest-neighbour
#' matching on the propensity with a caliper, and entropy balancing,
#' covariate-balancing weights for the effect on the treated. [dmlx()]
#' gives the cross-fitted ones (augmented inverse probability weighting,
#' targeted maximum likelihood, the interactive model).
#'
#' The propensity is the logistic regression of `d` on an intercept and
#' `x` (the learner of [lgtreg()]) unless `ps` gives it. `tgt` picks the
#' average effect (`"ate"`), the effect on the treated (`"att"`) or on the
#' controls (`"atc"`).
#' * `"ipw"`: the weighted mean of the treated less that of the controls,
#'   each normalised to sum to one (Hajek), with weights `1 / e` and `1 /
#'   (1 - e)` (ATE), `1` and `e / (1 - e)` (ATT), `(1 - e) / e` and `1`
#'   (ATC).
#' * `"gcomp"`: least squares of `y` on `x` in each arm; the mean
#'   difference of the two predictions over the target units.
#' * `"sub"`: `nsub` subclasses cut at the quantiles of the propensity
#'   among the target units; the subclasses' differences of means weighted
#'   by their target units.
#' * `"psm"`: greedy `k`:1 matching without replacement on the propensity,
#'   the treated (ATT) or the controls (ATC) in decreasing propensity
#'   ('MatchIt''s default order), within `calip` standard deviations of the
#'   propensity. The estimate is the mean over the matched units of the
#'   treated outcome less its controls' mean; `med` is their median, the
#'   median causal effect of the matched pairs; the standard error is that
#'   of a mean of independent matched sets.
#' * `"ebal"`: entropy balancing (Hainmueller 2012): control weights
#'   proportional to `exp(l'x)`, the nearest to uniform in entropy whose
#'   weighted means equal the treated means exactly; the ATT is the treated
#'   mean less the weighted control mean.
#'
#' Standard errors (but `"psm"`) are the standard deviations over `reps`
#' bootstrap samples of the units, with the propensity fitted again in
#' each; intervals are normal.
#'
#' A `clpre` list from [clprep()] may be given as `y`: its elements fill
#' the arguments by name (the `otis-*-cbatt` presets).
#'
#' @param y The outcome: finite numbers, one per unit.
#' @param d The treatment: 0 or 1 (or `FALSE` or `TRUE`), two units in each
#'   arm at least.
#' @param x The covariates: a numeric matrix, a data frame or list of
#'   numeric columns, or one numeric vector, one row per unit.
#' @param est `"ipw"`, `"gcomp"`, `"sub"`, `"psm"` or `"ebal"`.
#' @param tgt `"ate"`, `"att"` or `"atc"`.
#' @param ps `NULL` or the propensity of each unit, strictly between 0 and
#'   1.
#' @param k The controls matched to each unit: a whole number from 1.
#' @param calip `NULL` or the caliper in standard deviations of the
#'   propensity.
#' @param nsub The subclasses: a whole number from 2.
#' @param reps The bootstrap samples: a whole number from 0 (none).
#' @param seed The seed of the bootstrap: `NULL` to draw it from R's
#'   generator, or a whole number from 0.
#' @param level The confidence level of the normal interval.
#' @return An object of class `cbatt`, a list: `res`, one row (`est`, `se`,
#'   `lo`, `hi`, `med`, `n1`, `n0`); `ps`, the propensities (`NA` for
#'   `"gcomp"` and `"ebal"`); `w`, the weights (`"ipw"`, `"ebal"`; else
#'   `NULL`); `match` (`"psm"`: the rows of each `focal` unit and its
#'   `match`; else `NULL`); `est`; `tgt`; and `level`.
#' @references
#'   Hainmueller, J. (2012). Entropy balancing for causal effects.
#'   *Political Analysis*, 20(1), 25-46. \doi{10.1093/pan/mpr025}
#'
#'   Imbens, G. W. and Rubin, D. B. (2015). *Causal Inference for
#'   Statistics, Social, and Biomedical Sciences*. Cambridge University
#'   Press. \doi{10.1017/CBO9781139025751}
#'
#'   Austin, P. C. (2011). An introduction to propensity score methods for
#'   reducing the effects of confounding in observational studies.
#'   *Multivariate Behavioral Research*, 46(3), 399-424.
#'   \doi{10.1080/00273171.2011.568786}
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400), 200)
#' d <- rbinom(200, 1, plogis(x[, 1]))
#' y <- 2 * d + x[, 1] - x[, 2] + rnorm(200)
#' cbatt(y, d, x, reps = 99, seed = 1)$res
#' cbatt(y, d, x, est = "ebal", tgt = "att", reps = 99, seed = 1)$res
#' cbatt(y, d, x, est = "psm", tgt = "att", calip = 0.2)$res
#' @export
cbatt <- function(y, d, x, est = "ipw", tgt = "ate", ps = NULL, k = 1,
                  calip = NULL, nsub = 5, reps = 499, seed = NULL,
                  level = 0.95) {
  if (inherits(y, "clpre")) {
    return(.prep(cbatt, y, match.call(), parent.frame()))
  }
  .Call(C_clab_r_cbatt, y, d, x, est, tgt, ps, k, calip, nsub, reps, seed,
        level)
}

#' @export
print.cbatt <- function(x, ...) {
  cat("cbatt: ", x[["est"]], ", ", toupper(x[["tgt"]]), "\n", sep = "")
  print(x[["res"]], row.names = FALSE)
  invisible(x)
}

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

#' E-values: how strong an unmeasured confounder must be
#'
#' The E-value of VanderWeele and Ding (2017): the smallest strength of
#' association, on the risk-ratio scale, that an unmeasured confounder
#' would need with both the treatment and the outcome to explain away an
#' estimate; and that of the confidence limit nearer the null.
#'
#' For a risk ratio `RR >= 1` the E-value is `RR + sqrt(RR (RR - 1))`; a
#' ratio below 1 is inverted first. An odds or hazard ratio is a risk ratio
#' when the outcome is rare (`rare = TRUE`); otherwise an odds ratio is
#' taken to `sqrt(OR)` and a hazard ratio to `(1 - 0.5^sqrt(HR)) / (1 -
#' 0.5^sqrt(1 / HR))`. A standardised mean difference `d` is taken to
#' `exp(0.91 d)`, its limits to `exp(0.91 d -+ 1.78 sd)` from its standard
#' error `sd`. The limit's E-value is 1 when the interval holds the null.
#' These are the conversions of the 'EValue' package.
#'
#' @param est The estimate: one positive ratio, or one standardised
#'   difference.
#' @param lo,hi `NULL` or the confidence limits of a ratio.
#' @param type `"rr"`, `"or"`, `"hr"` or `"smd"`.
#' @param rare `TRUE` when the outcome is rare (an odds or hazard ratio is
#'   then read as a risk ratio).
#' @param sd `NULL` or the standard error of a standardised difference.
#' @return A named numeric vector: `rr` (the estimate as a risk ratio),
#'   `evalue` and `evlim` (the limit's E-value; `NA` without limits).
#' @references
#'   VanderWeele, T. J. and Ding, P. (2017). Sensitivity analysis in
#'   observational research: introducing the E-value. *Annals of Internal
#'   Medicine*, 167(4), 268-274. \doi{10.7326/M16-2607}
#' @examples
#' evlx(3.9, 1.8, 8.7)
#' evlx(0.8, 0.71, 0.91, type = "hr")
#' evlx(0.5, type = "smd", sd = 0.1)
#' @export
evlx <- function(est, lo = NULL, hi = NULL, type = "rr", rare = FALSE,
                 sd = NULL) {
  .Call(C_clab_r_evlx, est, lo, hi, type, rare, sd)
}
