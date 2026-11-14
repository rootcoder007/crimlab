#' Calibration line of an assay
#'
#' Fits the calibration line of a quantitative assay (the response, an
#' analyte to internal standard peak-area ratio say, against the known
#' concentration of each calibrator) by weighted least squares, and reads
#' unknown samples off it. Weights `1/x^2` or `1/x` keep the precise low
#' calibrators from being swamped by the high ones, as is usual when the
#' error grows with the concentration. The limits of detection and of
#' quantitation follow the calibration-curve rule of ICH Q2: `k[1]` and
#' `k[2]` (3.3 and 10) times the residual standard error over the slope.
#'
#' @param conc The concentrations of the calibrators: positive numbers, at
#'   least three.
#' @param resp The response of each calibrator: finite numbers.
#' @param w The weights: `"1/x^2"`, `"1/x"`, `"none"`, or positive numbers,
#'   one per calibrator.
#' @param unk `NULL`, or the responses of unknown samples: numbers or `NA`.
#' @param k The multipliers of the limits of detection and quantitation:
#'   two positive numbers.
#' @return A list: `fit`, a named vector (`int` and `slope` of the line,
#'   `r2` the weighted R squared, `sigma` the residual standard error on
#'   `n - 2` degrees of freedom, `lod` and `loq` the limits); `unk`, `NULL`
#'   or a data frame with one row per unknown (`resp`, `conc` the
#'   concentration read off the line, `flag`: `below_lod`, `below_loq` or
#'   `quantifiable`); `w`, the weighting.
#' @references
#'   ICH (2023). *ICH Harmonised Guideline Q2(R2): Validation of
#'   Analytical Procedures*, final version adopted 1 November 2023, section
#'   3.2.3.3. International Council for Harmonisation of Technical
#'   Requirements for Pharmaceuticals for Human Use.
#'   <https://www.ich.org/page/quality-guidelines>
#' @examples
#' conc <- c(5, 10, 25, 50, 100, 250)
#' resp <- c(0.051, 0.098, 0.26, 0.49, 1.02, 2.47)
#' txcal(conc, resp, unk = c(0.03, 0.6))
#' @export
txcal <- function(conc, resp, w = "1/x^2", unk = NULL, k = c(3.3, 10)) {
  .Call(C_clab_r_txcal, conc, resp, w, unk, k)
}

#' Values below the limit of detection
#'
#' Substitutes the values of an assay that fall below the limit of
#' detection (left-censored) by a fraction of that limit: `lod / 2`,
#' `lod / sqrt(2)` (Hornung and Reed 1990) or `lod` itself. A value not
#' measured (`NA`) stays `NA` and is not counted as censored. Substitution
#' is the simple choice; with much censoring it biases the mean and the
#' spread.
#'
#' @param x The measured values: numbers, `NA` for not measured. A `clpre`
#'   list from [clprep()] is taken too (the `tox-lodimp` preset): its
#'   elements fill the arguments by name.
#' @param lod The limit of detection: one positive number, or one per value
#'   (each sample's own limit).
#' @param meth `"half"`, `"sqrt2"` or `"lod"`.
#' @return A list: `tab`, a data frame with one row per value (`x`, `lod`,
#'   `imp` the value after substitution, `cens` whether it was below the
#'   limit); `frac`, the share of the measured values below the limit
#'   (`NA` when none was measured); `nmiss`, the values not measured.
#' @references
#'   Hornung, R. W. and Reed, L. D. (1990). Estimation of average
#'   concentration in the presence of nondetectable values. *Applied
#'   Occupational and Environmental Hygiene*, 5(1), 46-51.
#'   \doi{10.1080/1047322X.1990.10389587}
#' @examples
#' lodimp(c(0.4, 0.05, NA, 1.2, 0.08), lod = 0.1)
#' @export
lodimp <- function(x, lod, meth = "half") {
  if (inherits(x, "clpre")) {
    return(.prep(lodimp, x, match.call(), parent.frame()))
  }
  .Call(C_clab_r_lodimp, x, lod, meth)
}

#' Postmortem redistribution ratio
#'
#' After death drugs move from organs into the central blood vessels, so
#' that heart blood can read well above the concentration at death. The
#' ratio of the central (heart) to the peripheral (femoral) concentration
#' flags it: at most `cut[1]` (1) minimal, at most `cut[2]` (2) modest,
#' above that significant, when the peripheral value is the better guide
#' to the concentration at death.
#'
#' @param cen,per The central and peripheral blood concentrations: positive
#'   numbers, one per case.
#' @param cut The two cut points of the ratio: positive, increasing.
#' @return A data frame with one row per case: `cen`, `per`, `ratio` (`cen /
#'   per`) and `redist`, a factor (`minimal`, `modest`, `significant`).
#' @references
#'   Yarema, M. C. and Becker, C. E. (2005). Key concepts in postmortem
#'   drug redistribution. *Clinical Toxicology*, 43(4), 235-241.
#'   \doi{10.1081/CLT-58950}
#' @examples
#' pmrrat(c(0.4, 1.5, 0.9), c(0.5, 1, 0.2))
#' @export
pmrrat <- function(cen, per, cut = c(1, 2)) {
  .Call(C_clab_r_pmrrat, cen, per, cut)
}

#' Ethanol in a decomposed body: drunk before death or made after
#'
#' Microbes make ethanol in a body after death, and with it the higher
#' alcohols n-propanol and n-butanol; ethyl glucuronide, a metabolite, is
#' made only by a living body. So for each case: ethyl glucuronide present
#' (above 0) says the ethanol was taken before death; otherwise a higher
#' alcohol above `thr` says it was (at least partly) made after death;
#' otherwise the signals do not decide.
#'
#' @param eth The ethanol concentration: numbers from 0, one per case.
#' @param prop,but The n-propanol and n-butanol concentrations: numbers or
#'   `NA` (not measured), one per case or one for all.
#' @param etg The ethyl glucuronide concentration: numbers or `NA` (not
#'   measured), one per case or one for all.
#' @param thr The level a higher alcohol must exceed: one number from 0.
#' @return A data frame with one row per case: `eth`, `prop`, `but`, `etg`,
#'   `higher` (the larger higher alcohol, `NA` when neither was measured)
#'   and `verdict`, a factor (`antemortem`, `postmortem`,
#'   `indeterminate`).
#' @references
#'   Kugelberg, F. C. and Jones, A. W. (2007). Interpreting results of
#'   ethanol analysis in postmortem specimens: a review of the literature.
#'   *Forensic Science International*, 165(1), 10-29.
#'   \doi{10.1016/j.forsciint.2006.05.004}
#' @examples
#' etohlr(c(0.8, 0.5, 0.3), prop = c(0, 0.02, NA), etg = c(1.1, NA, NA))
#' @export
etohlr <- function(eth, prop = 0, but = 0, etg = NA, thr = 0) {
  .Call(C_clab_r_etohlr, eth, prop, but, etg, thr)
}

#' Reliability of specimen matrices in a submerged or decomposed body
#'
#' Ranks the specimen matrices of a case (vitreous humour, liver, femoral
#' blood) by how far their results can be trusted. Each matrix starts from
#' a reliability in a fresh body, `base`, which the user sets from the
#' literature for the analyte at hand; decomposition lowers every matrix by
#' `stage / smax`, at most by `scap`, and water dilutes the matrices marked
#' `dil` (blood, urine) by a further `days / dmax`, at most `dcap`. A
#' transparent ordering, not an estimate.
#'
#' @param base The reliability of each matrix in a fresh body: numbers from
#'   0, named by matrix.
#' @param dil Whether water dilutes each matrix: `TRUE` or `FALSE`, one per
#'   matrix.
#' @param days The days the body was submerged: one number from 0.
#' @param stage The decomposition stage (0 fresh): one number from 0.
#' @param dmax,smax The days and the stage at which the loss reaches its
#'   rate's full scale: positive numbers.
#' @param dcap,scap The largest losses to dilution and to decomposition:
#'   numbers from 0 to 1.
#' @return A data frame with one row per matrix, the most reliable first:
#'   `matrix` (the names of `base`, or their numbers), `base`, `rel` the
#'   reliability after the losses, `dil`, `rank`.
#' @references
#'   Dinis-Oliveira, R. J., Carvalho, F., Duarte, J. A., Remiao, F.,
#'   Marques, A., Santos, A. and Magalhaes, T. (2010). Collection of
#'   biological samples in forensic toxicology. *Toxicology Mechanisms and
#'   Methods*, 20(7), 363-414. \doi{10.3109/15376516.2010.497976}
#' @examples
#' base <- c(vitreous = 0.95, liver = 0.8, femoral = 0.55, heart = 0.35)
#' mtxrel(base, dil = c(FALSE, FALSE, TRUE, TRUE), days = 12, stage = 3)
#' @export
mtxrel <- function(base, dil, days = 0, stage = 0, dmax = 30, dcap = 0.8,
                   smax = 10, scap = 0.6) {
  .Call(C_clab_r_mtxrel, base, dil, days, stage, dmax, dcap, smax, scap)
}
