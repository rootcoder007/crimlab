#' Machine-checked proofs behind the package
#'
#' Names the 'Lean' theorem each function rests on. Theorems with source
#' `"rmorie"` are part of the 'rmorie' research ledger, and `file` is the path
#' inside its `research/lean` directory. Theorems with source `"csls"` ship
#' with this package, and `file` is the path below
#' `system.file(package = "csls")`.
#'
#' @param fn `NULL` for every function, or function names.
#' @return A data frame: `fn`, `theorem`, `src` and `file`.
#' @examples
#' proof("vodtst")
#' @export
proof <- function(fn = NULL) {
  row <- function(fn, ns, thm, src, file) {
    data.frame(fn = fn, theorem = paste0(ns, ".", thm), src = src,
               file = file)
  }
  led <- function(fn, ns, thm, file) {
    row(fn, ns, thm, "rmorie", paste0("Researchproofs/", file))
  }
  own <- function(fn, thm) {
    row(fn, "Csls.Forensics", thm, "csls", "lean/Forensics.lean")
  }
  t <- rbind(
    led("vodtst", "Research.P2.Veil",
        c("odds_ratio_identifies_bias", "share_sign",
          "exposure_shift_confounds"), "P2Veil.lean"),
    led("orrisk", "Research.P2",
        c("or_eq_rr_mul", "rr_between", "or_overstates"),
        "P2RelativeRisk.lean"),
    led("outtst", "Research.P2.OutcomeTest",
        c("infra_marginality", "equal_thresholds_unequal_hits"),
        "P2OutcomeTest.lean"),
    led("dspbm", "Research.P2",
        c("benchmark_product", "benchmark_not_additive"), "P2Benchmark.lean"),
    led("dyadrt", "Research.P2",
        c("null_slope_positive", "pair_exposure_rate_constant"),
        "P2Interracial.lean"),
    own("gtmp", c("bn_hom_zero", "bn_het_zero", "bn_hom_ge_product",
                  "bn_het_not_ge_product")),
    own("lrpost", c("post_eq", "post_lt_half")),
    own("cllr", "cllr_neutral"),
    own("dbsrch", "np_rule"),
    own("mixlr", "mix_two_alleles"),
    own("kinidx", "po_identical_het"),
    own("patidx", "post_even_prior")
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
