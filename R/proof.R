#' Machine-checked proofs behind the package
#'
#' Names the 'Lean' theorem each function rests on. The theorems ship with
#' the package: `file` is the path below `system.file(package = "crimlab")`, and
#' each file compiles against 'Mathlib'.
#'
#' @param fn `NULL` for every function, or function names.
#' @return A data frame: `fn`, `theorem` and `file`.
#' @examples
#' proof("vodtst")
#' @export
proof <- function(fn = NULL) {
  .Call(C_clab_r_proof, fn)
}
