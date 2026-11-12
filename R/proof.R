#' Machine-checked proofs behind the package
#'
#' Names the 'Lean' theorem each function rests on. The theorems are part of
#' the 'rmorie' research ledger; the file is the path inside its
#' `research/lean` directory.
#'
#' @param fn `NULL` for every function, or function names.
#' @return A data frame: `fn`, `theorem` and `file`.
#' @examples
#' proof("vodtst")
#' @export
proof <- function(fn = NULL) {
  p <- "Researchproofs/"
  t <- data.frame(
    fn = c("vodtst", "vodtst", "vodtst", "orrisk", "orrisk", "orrisk", "outtst",
           "outtst", "dspbm", "dspbm", "dyadrt", "dyadrt"),
    theorem = c(
      "Research.P2.Veil.odds_ratio_identifies_bias",
      "Research.P2.Veil.share_sign",
      "Research.P2.Veil.exposure_shift_confounds",
      "Research.P2.or_eq_rr_mul",
      "Research.P2.rr_between",
      "Research.P2.or_overstates",
      "Research.P2.OutcomeTest.infra_marginality",
      "Research.P2.OutcomeTest.equal_thresholds_unequal_hits",
      "Research.P2.benchmark_product",
      "Research.P2.benchmark_not_additive",
      "Research.P2.null_slope_positive",
      "Research.P2.pair_exposure_rate_constant"
    ),
    file = paste0(p, c("P2Veil.lean", "P2Veil.lean", "P2Veil.lean",
                       "P2RelativeRisk.lean", "P2RelativeRisk.lean",
                       "P2RelativeRisk.lean",
                       "P2OutcomeTest.lean", "P2OutcomeTest.lean",
                       "P2Benchmark.lean", "P2Benchmark.lean",
                       "P2Interracial.lean", "P2Interracial.lean"))
  )
  if (is.null(fn)) return(t)
  if (!is.character(fn) || anyNA(fn)) {
    stop("`fn` must be function names.", call. = FALSE)
  }
  miss <- setdiff(fn, t[["fn"]])
  if (length(miss)) {
    stop("no proof recorded for: ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  out <- t[t[["fn"]] %in% fn, , drop = FALSE]
  rownames(out) <- NULL
  out
}
