#' Machine-checked proofs behind the package
#'
#' Names the 'Lean' theorem each function rests on. The theorems ship with
#' the package: `file` is the path below `system.file(package = "csls")`, and
#' each file compiles against 'Mathlib'.
#'
#' @param fn `NULL` for every function, or function names.
#' @return A data frame: `fn`, `theorem` and `file`.
#' @examples
#' proof("vodtst")
#' @export
proof <- function(fn = NULL) {
  row <- function(fn, mod, thm) {
    data.frame(fn = fn, theorem = paste0("Csls.", mod, ".", thm),
               file = paste0("lean/", mod, ".lean"))
  }
  t <- rbind(
    row("vodtst", "Veil", c("odds_ratio_identifies_bias", "share_sign",
                            "exposure_shift_confounds")),
    row("orrisk", "RelRisk", c("or_eq_rr_mul", "rr_between", "or_overstates")),
    row("outtst", "OutcomeTest", c("infra_marginality",
                                   "equal_thresholds_unequal_hits")),
    row("dspbm", "Benchmark", c("benchmark_product",
                                "benchmark_not_additive")),
    row("dyadrt", "Dyad", c("null_slope_positive",
                            "pair_exposure_rate_constant")),
    row("gtmp", "Forensics", c("bn_hom_zero", "bn_het_zero",
                               "bn_hom_ge_product", "bn_het_not_ge_product")),
    row("lrpost", "Forensics", c("post_eq", "post_lt_half")),
    row("cllr", "Forensics", "cllr_neutral"),
    row("dbsrch", "Forensics", "np_rule"),
    row("mixlr", "Forensics", "mix_two_alleles"),
    row("kinidx", "Forensics", "po_identical_het"),
    row("patidx", "Forensics", "post_even_prior"),
    row("clrsmp", "Clearance", c("agg_eq_weighted", "kitagawa",
                                 "change_ge_mix", "fixed_mix_no_reversal",
                                 "simpson_witness")),
    row("expbnd", "Exposure", c("rate_bounds", "rate_lower_attained",
                                "rate_upper_attained", "disparity_bounds",
                                "disparity_sign_identified")),
    row("misrd", "Misclass", c("binary_attenuation", "binary_sign")),
    row("rtmhs", "RegMean", c("selected_change_nonpos",
                              "low_selected_change_nonneg"))
  )
  if (is.null(fn)) return(t)
  if (!is.character(fn) || anyNA(fn)) {
    stop("`fn` must be function names.", call. = FALSE)
  }
  miss <- setdiff(fn, t[["fn"]])
  if (length(miss)) {
    stop("no proof recorded for: ", toString(miss), ".",
         call. = FALSE)
  }
  out <- t[t[["fn"]] %in% fn, , drop = FALSE]
  rownames(out) <- NULL
  out
}
