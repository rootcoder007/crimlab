#' Eyewitness lineups: ROC curves, partial areas and confidence-accuracy
#'
#' Asks how well a lineup procedure lets witnesses tell the culprit from an
#' innocent suspect, and how much a witness's confidence is worth. Each
#' witness saw a lineup that held the culprit (target present, `tp`) or did
#' not (target absent), picked the suspect, a filler or nobody (`pick`),
#' and stated a confidence (`conf`, higher is surer). Counting only the
#' identifications made at confidence `c` or above gives, at each cut `c`
#' from the highest down, the hit rate (suspect identifications in
#' target-present lineups over those lineups) and the false alarm rate
#' (innocent-suspect identifications in target-absent lineups over those
#' lineups). The points trace the receiver operating characteristic (ROC)
#' of Mickes, Flowe and Wixted (2012); a procedure whose curve lies higher
#' discriminates better. When target-absent lineups have no designated
#' innocent suspect, give the lineup `size`: every identification in a
#' target-absent lineup then counts, divided by `size` (the rate at which
#' an innocent suspect would be picked in a fair lineup).
#'
#' * `roc`: the hit rate `hr`, false alarm rate `far` and diagnosticity
#'   ratio `dr = hr / far` at each cut.
#' * `pauc`: the partial area under the curve from a false alarm rate of 0
#'   to `far` (by default the smallest of the procedures' largest false
#'   alarm rates, so all are compared over the same range), by the
#'   trapezoid rule; a curve that ends before `far` is taken flat from its
#'   last point. With `nboot` the standard errors come from resampling the
#'   target-present and target-absent witnesses of each procedure, and
#'   each procedure is compared with the first by `z = diff / dse`, as
#'   `pROC::roc.test(method = "bootstrap")` compares partial areas.
#' * `cac`: the confidence-accuracy characteristic (Mickes 2015): at each
#'   confidence level alone, the share of suspect identifications that are
#'   correct, `hr / (hr + far)` of that level, which assumes guilty and
#'   innocent suspects equally common.
#'
#' @param tp Whether each witness's lineup held the culprit: 0 or 1 (or
#'   `FALSE` or `TRUE`).
#' @param pick What each witness picked: `"suspect"`, `"filler"` or
#'   `"none"`.
#' @param conf The confidence of each witness: numbers, higher for surer;
#'   `NA` is allowed where `pick` is `"none"`.
#' @param proc `NULL` or the procedure of each witness (simultaneous,
#'   sequential): labels.
#' @param size `NULL` when target-absent lineups hold a designated innocent
#'   suspect, or the number of people in a lineup: one whole number from 2.
#' @param far `NULL` or the false alarm rate up to which the partial area
#'   runs: one number above 0, at most 1.
#' @param nboot The bootstrap samples: one whole number from 0 (none).
#' @param seed The seed of the bootstrap: `NULL` to draw it from R's
#'   generator, or one whole number from 0.
#' @param level The confidence level of the intervals: one number strictly
#'   between 0 and 1.
#' @return An object of class `lnroc`, a list: `roc`, one row per
#'   procedure (`proc`, when given) and cut (`conf`, `hr`, `far`, `dr`);
#'   `cac`, one row per procedure and confidence level (`conf`, `hr`,
#'   `far`, `acc`); `pauc`, one row per procedure (`pauc`, `se`, `lo`,
#'   `hi`, and `diff`, `dse`, `p` against the first procedure; the
#'   standard errors `NA` without the bootstrap); `far`, the false alarm
#'   rate used; `size`; and `nboot`.
#' @references
#'   Mickes, L., Flowe, H. D. and Wixted, J. T. (2012). Receiver operating
#'   characteristic analysis of eyewitness memory: comparing the diagnostic
#'   accuracy of simultaneous versus sequential lineups. *Journal of
#'   Experimental Psychology: Applied*, 18(4), 361-376.
#'   \doi{10.1037/a0030609}
#'
#'   Mickes, L. (2015). Receiver operating characteristic analysis and
#'   confidence-accuracy characteristic analysis in investigations of system
#'   variables and estimator variables that affect eyewitness memory.
#'   *Journal of Applied Research in Memory and Cognition*, 4(2), 93-102.
#'   \doi{10.1016/j.jarmac.2015.01.003}
#' @seealso [lnsize()] for the fairness of a lineup.
#' @examples
#' set.seed(3)
#' n <- 400
#' tp <- rep(0:1, each = n / 2)
#' proc <- rep(c("sim", "seq"), n / 2)
#' ch <- c("suspect", "filler", "none")
#' pick <- ifelse(tp == 1, sample(ch, n, TRUE, c(0.5, 0.15, 0.35)),
#'                sample(ch, n, TRUE, c(0.1, 0.3, 0.6)))
#' conf <- ifelse(pick == "none", NA, sample(1:5, n, TRUE))
#' r <- lnroc(tp, pick, conf, proc = proc, nboot = 200, seed = 1)
#' r$pauc
#' plot(r$roc$far, r$roc$hr, col = factor(r$roc$proc))
#' @export
lnroc <- function(
  tp,
  pick,
  conf,
  proc = NULL,
  size = NULL,
  far = NULL,
  nboot = 0,
  seed = NULL,
  level = 0.95
) {
  .Call(C_clab_r_lnroc, tp, pick, conf, proc, size, far, nboot, seed, level)
}

#' @export
print.lnroc <- function(x, ...) {
  cat(
    "lnroc: partial area up to a false alarm rate of ",
    format(signif(x[["far"]], 4)),
    "\n",
    sep = ""
  )
  print(x[["pauc"]], row.names = FALSE)
  cat("ROC points:\n")
  print(x[["roc"]], row.names = FALSE)
  invisible(x)
}

#' Fairness of a lineup from mock witnesses
#'
#' Asks whether a lineup singles out its suspect. People who never saw the
#' crime (mock witnesses) read the witness's description and pick the
#' lineup member who fits it best; in a fair lineup of `k` they pick each
#' member about equally often. From their choices:
#'
#' * `psusp`, the share who picked the suspect, and `pbin`, the binomial
#'   chance of that many or more if each member were picked with chance
#'   `1 / k` (a small `pbin` is bias against the suspect);
#' * `fsize`, the functional size, the number of mock witnesses over those
#'   who picked the suspect (Wells, Leippe and Ostrom 1979; `k` for a fair
#'   lineup);
#' * `esize`, Malpass's (1981) effective size, `k - sum |o_i - e| / (2 e)`
#'   with `o_i` the choices of member `i` and `e = n / k`: how many members
#'   are plausible;
#' * `etr`, Tredoux's (1998) effective size `1 / sum p_i^2`, `p_i` the share
#'   choosing member `i`, with a delta-method interval from the multinomial
#'   variance of `sum p_i^2`, `4 (sum p^3 - (sum p^2)^2) / n`.
#'
#' @param choice The lineup position each mock witness picked: whole
#'   numbers from 1.
#' @param k `NULL` for the largest position picked, or the lineup's size:
#'   one whole number from 2.
#' @param susp The suspect's position: one whole number from 1 to `k`.
#' @inheritParams lnroc
#' @return An object of class `lnsize`, a list: `count`, the choices of each
#'   position; `stat`, a named vector (`n`, `k`, `psusp`, `pbin`, `fsize`,
#'   `esize`, `etr`, `etrlo`, `etrhi`); and `susp`.
#' @references
#'   Malpass, R. S. (1981). Effective size and defendant bias in eyewitness
#'   identification lineups. *Law and Human Behavior*, 5(4), 299-309.
#'   \doi{10.1007/BF01044945}
#'
#'   Tredoux, C. G. (1998). Statistical inference on measures of lineup
#'   fairness. *Law and Human Behavior*, 22(2), 217-237.
#'   \doi{10.1023/A:1025746220886}
#'
#'   Wells, G. L., Leippe, M. R. and Ostrom, T. M. (1979). Guidelines for
#'   empirically assessing the fairness of a lineup. *Law and Human
#'   Behavior*, 3(4), 285-293. \doi{10.1007/BF01039807}
#' @examples
#' # 40 mock witnesses, a lineup of six, the suspect in position 3
#' ch <- c(rep(3, 15), rep(1, 5), rep(2, 6), rep(4, 4), rep(5, 7), rep(6, 3))
#' lnsize(ch, k = 6, susp = 3)
#' @export
lnsize <- function(choice, k = NULL, susp = 1, level = 0.95) {
  .Call(C_clab_r_lnsize, choice, k, susp, level)
}

#' @export
print.lnsize <- function(x, ...) {
  s <- x[["stat"]]
  cat(
    "lnsize: ",
    s[["n"]],
    " mock witnesses, lineup of ",
    s[["k"]],
    ", suspect in position ",
    x[["susp"]],
    "\n",
    sep = ""
  )
  print(s)
  invisible(x)
}
