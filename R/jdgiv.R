#' The judge-leniency instrument
#'
#' When cases are assigned to judges at random, the judge's leniency is an
#' instrument for the decision (detention, a custodial sentence): the
#' difference in outcomes between lenient and strict judges over the
#' difference in decisions is the Wald ratio (Imbens and Angrist 1994).
#' It is the average effect on the compliers, the cases whose decision
#' depends on the judge, only if no case defies (is detained by the
#' lenient judge and released by the strict one). With a share of defiers
#' `s` whose effect is `e`, the compliers number the first stage plus `s`
#' and their effect is `(wald * fs + s * e) / (fs + s)`; `sens` gives it
#' over a grid, so a reader sees how far the conclusion rests on
#' monotonicity (de Chaisemartin 2017). Code `z` so that the first stage
#' is positive (1 for the judge whose decision is more often 1); with a
#' negative one the complier share is negative and `ceff` is `NA`.
#'
#' Given the decisions and outcomes of a whole population under both
#' judges (`z = NULL`, `d` and `y` two columns each), the exact
#' decomposition: the shares and mean effects of compliers, defiers,
#' always- and never-takers, the intention-to-treat effect and the first
#' stage, whose ratio is the compliers' effect when there are no defiers.
#' See `proof("jdgiv")`.
#'
#' @param z `NULL`, or the instrument of each case: 0 for the strict judge,
#'   1 for the lenient (or `FALSE`, `TRUE`).
#' @param d The decision of each case: 0 or 1; with `z = NULL`, two
#'   columns (a matrix or data frame), the decisions under the strict and
#'   under the lenient judge.
#' @param y The outcome of each case: numbers; with `z = NULL`, two
#'   columns, the outcomes without and with the decision.
#' @param w `NULL`, or a weight from 0 for each case.
#' @param dsh Defier shares for the sensitivity grid: shares from 0 to 1.
#' @param deff Defier effects for the sensitivity grid: numbers.
#' @return An object of class `jdgiv`, a list: `type` (`"obs"` or
#'   `"pop"`); `est`, a named vector (observed: `fs` the first stage, `itt`,
#'   `wald`; population: `itt`, `fs`, `wald`, `ate`); `shares` (observed:
#'   `complier`, `always`, `never` if no one defies; population: also
#'   `defier`); and `sens` (observed: `dsh`, `deff`, `cshare` the complier
#'   share, `ceff` their effect) or `effects` (population: the mean effect
#'   of each type, `NA` for an empty one).
#' @references
#'   Imbens, G. W. and Angrist, J. D. (1994). Identification and
#'   estimation of local average treatment effects. *Econometrica*, 62(2),
#'   467-475. \doi{10.2307/2951620}
#'
#'   de Chaisemartin, C. (2017). Tolerating defiance? Local average
#'   treatment effects without monotonicity. *Quantitative Economics*,
#'   8(2), 367-396. \doi{10.3982/QE601}
#' @examples
#' set.seed(7)
#' z <- rbinom(1000, 1, 0.5)
#' d <- rbinom(1000, 1, ifelse(z == 1, 0.3, 0.5))
#' y <- rbinom(1000, 1, 0.2 + 0.1 * d)
#' jdgiv(z, d, y)
#' @export
jdgiv <- function(z = NULL, d, y, w = NULL, dsh = c(0, 0.05, 0.1), deff = 0) {
  structure(.Call(C_clab_r_jdgiv, z, d, y, w, dsh, deff), class = "jdgiv")
}

#' Do judges' leniency rankings hold? Monotonicity checks
#'
#' The decision rate of each judge (with a chi-square test that the rate
#' is the same for all, and its range), and two checks of the
#' monotonicity a judge-leniency design rests on.
#'
#' * With outcomes `y` in `[lo, hi]`: if judges are nested (a stricter
#'   judge detains everyone a more lenient one detains), the outcome means
#'   of two judges differ by at most `(hi - lo)` times the difference in
#'   their rates. Every pair is checked; a violation refutes the nesting
#'   of that pair. The ratio `late` is the effect on the marginal cases
#'   between them.
#' * With groups of cases `grp`: nested judges are ordered the same way in
#'   every group. For each pair of judges with at least `minn` cases in a
#'   group, the stricter judge overall should not be more lenient there;
#'   `z` is the pooled two-proportion statistic of the lenient judge's
#'   group rate over the stricter's, its one-sided p-value is adjusted
#'   over all comparisons, and `slopes` regresses each group's rates on the
#'   overall ones (Frandsen, Lefgren and Leslie 2023).
#'
#' See `proof("jmono")`.
#'
#' @param judge The judge of each case: labels.
#' @param d The decision of each case: 0 or 1 (or `FALSE`, `TRUE`).
#' @param y `NULL`, or the outcome of each case: numbers.
#' @param grp `NULL`, or the group of each case: labels.
#' @param w `NULL`, or a weight from 0 for each case (for the rates and
#'   outcome means; the tests use the counts).
#' @param lo,hi `NULL` (the smallest and largest outcome) or the bounds of
#'   the outcome.
#' @param minn The fewest cases a judge needs in a group to be compared
#'   there: one number from 0.
#' @param adj The adjustment of the p-values: `"none"`, `"bonferroni"`,
#'   `"holm"`, `"hochberg"`, `"BH"` or `"BY"`.
#' @param alpha The level at which an adjusted p-value marks a reversal as
#'   significant.
#' @return An object of class `jmono`, a list: `judges` (`judge`, `n`,
#'   `rate`, with `y` `outcome`; most lenient first); `homog` (`chisq`,
#'   `df`, `p`, `min`, `max`); `slope`, with `y`, one row per pair (`j`
#'   the more lenient, `k`, `dp`, `dy`, `bound`, `viol`, `late`); `sub`,
#'   with `grp`, a list of `pairs` (`group`, `stricter`, `lenient`, their
#'   overall rates `rates`, `ratel`, their group rates `grates`, `gratel`,
#'   their cases `ns`, `nl`, `z`, `p`, `padj`, `reversed`, `signif`),
#'   `summary` (`pairs`, `reversed`, `share`, `signif`, `refuted` the pairs
#'   significantly reversed, `judges` in them) and `slopes` (`group`,
#'   `judges`, `cases`, `slope`); and `alpha`. The pairs grow as the square
#'   of the judges.
#' @references Frandsen, B. R., Lefgren, L. J. and Leslie, E. C. (2023).
#'   Judging judge fixed effects. *American Economic Review*, 113(1),
#'   253-277. \doi{10.1257/aer.20201860}
#' @examples
#' set.seed(9)
#' j <- sample(paste0("J", 1:6), 1200, TRUE)
#' len <- c(J1 = 0.2, J2 = 0.3, J3 = 0.35, J4 = 0.4, J5 = 0.5, J6 = 0.6)
#' d <- rbinom(1200, 1, len[j])
#' g <- sample(c("adult", "youth"), 1200, TRUE)
#' r <- jmono(j, d, grp = g)
#' r$sub$summary
#' @export
jmono <- function(
  judge,
  d,
  y = NULL,
  grp = NULL,
  w = NULL,
  lo = NULL,
  hi = NULL,
  minn = 30,
  adj = "holm",
  alpha = 0.05
) {
  r <- .Call(C_clab_r_jmono, judge, d, y, grp, w, lo, hi, minn, adj, alpha)
  structure(r, class = "jmono")
}

#' Sentences at, below and above a mandatory minimum
#'
#' The shares of sentences at the mandatory minimum, below it (where a
#' court found a way around it) and above it, and the mean excess of those
#' above, overall and by group (offence, court).
#'
#' @param sent The sentence of each case: numbers.
#' @param min The mandatory minimum: one number, or one per case.
#' @param grp `NULL`, or the group of each case: labels.
#' @return A data frame of class `mmshr` with one row per group and a last
#'   row `"all"`: `group`, `n`, `at`, `below`, `above` (shares) and
#'   `mexcess` (0 when none is above).
#' @examples
#' mmshr(c(365, 400, 300, 365, 730), 365, c("a", "a", "b", "b", "b"))
#' @export
mmshr <- function(sent, min, grp = NULL) {
  r <- .Call(C_clab_r_mmshr, sent, min, grp)
  structure(r, class = c("mmshr", "data.frame"))
}

#' @export
print.jdgiv <- function(x, ...) {
  cat("jdgiv:", if (x[["type"]] == "obs") "observed" else "population", "\n")
  print(x[["est"]])
  print(x[["shares"]])
  if (!is.null(x[["sens"]])) {
    print(x[["sens"]], row.names = FALSE)
  }
  if (!is.null(x[["effects"]])) {
    print(x[["effects"]])
  }
  invisible(x)
}

#' @export
print.jmono <- function(x, ...) {
  h <- x[["homog"]]
  cat(sprintf(
    "jmono: %d judges, rates %.3g to %.3g (chi-square p = %.3g)\n",
    nrow(x[["judges"]]),
    h[["min"]],
    h[["max"]],
    h[["p"]]
  ))
  if (!is.null(x[["slope"]])) {
    cat(sprintf(
      "slope test: %d of %d pairs violate the bound\n",
      sum(x[["slope"]][["viol"]]),
      nrow(x[["slope"]])
    ))
  }
  if (!is.null(x[["sub"]])) {
    print(x[["sub"]][["summary"]])
  }
  invisible(x)
}

#' @export
print.mmshr <- function(x, ...) {
  cat("mmshr: sentences against the minimum\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}
