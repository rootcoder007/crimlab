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
#'   trapezoid rule. A curve that ends before `far` is continued by the
#'   straight line from its last point to the point where every lineup
#'   counts as an identification (false alarm rate 1, or `1 / size` with
#'   `size`, and hit rate 1), as the package pROC that Mickes, Flowe and
#'   Wixted (2012) used joins a curve to (1, 1). With `nboot` the standard
#'   errors come from resampling the target-present and target-absent
#'   witnesses of each procedure, and each procedure is compared with the
#'   first by `z = diff / dse`, their `D` statistic, as
#'   `pROC::roc.test(method = "bootstrap")` compares partial areas.
#'   Their partial areas (0.13 and 0.09 in Experiment 1a, p. 371) are on
#'   the scale of all identifications in target-absent lineups, not
#'   divided by the lineup's size as the false alarm rates of their Table
#'   3 are: with `size`, `size * pauc` is on their scale.
#' * `cac`: the confidence-accuracy characteristic (Mickes 2015): at each
#'   confidence level alone, the share of suspect identifications that are
#'   correct, `hr / (hr + far)` of that level, which assumes guilty and
#'   innocent suspects equally common.
#'
#' The ROC points and the partial areas reproduce the rates of Table 3 of
#' Mickes, Flowe and Wixted (2012) and, on their scale, its partial areas
#' (see the tests).
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
#' @inheritParams cntreg
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
#' member about equally often. A choice of 0 records a mock witness who
#' chose no one. With `n` mock witnesses, `m` of whom chose a member, and
#' `o_i` the choices of member `i`:
#'
#' * `psusp`, the share of the `n` who picked the suspect, and `pbin`, the
#'   binomial chance of that many or more among the `m` if each member were
#'   picked with chance `1 / k` (a small `pbin` is bias against the
#'   suspect);
#' * `fsize`, the functional size of Wells, Leippe and Ostrom (1979), `n`
#'   over the number who picked the suspect (`k` for a fair lineup);
#' * `esize`, the effective size of Malpass (1981), `ka - sum |o_i - ea| /
#'   (2 ea)`: the members nobody chose (null foils) are left out, `ka` is
#'   the number of members chosen at least once, `ea = m / ka`, and the sum
#'   runs over those `ka` members; it reduces the size by how far members
#'   are chosen below chance. `esizeb` is the same without the adjustment
#'   (`k` and `e = m / k`, the sum over all members), the form Tredoux
#'   (1996, Chapter 6, equation 5; 1999, p. S13) prefers because a member
#'   can go unchosen by chance;
#' * `etr`, the effective size `E' = 1 / (1 - I) = 1 / sum p_i^2` of
#'   Tredoux (1998), `p_i = o_i / m`, with its interval: the large-sample
#'   interval of `I = 1 - sum p_i^2`, whose variance is `4 (sum p^3 - (sum
#'   p^2)^2) / m` (Agresti and Agresti 1978), with both ends taken to `1 /
#'   (1 - I)` (Tredoux 1996, Chapter 6, equations 8, 9 and 11). `etrhi` is
#'   `Inf` when the upper end of `I` reaches 1.
#'
#' The tests reproduce Tables 4.1 (p. 67), 4.2 (p. 70) and 6.6 (p. 144) of
#' Tredoux (1996) and Table 1 (p. S10) of Tredoux (1999). Two printed
#' values there do not follow from their own counts and formulas: lineup i
#' of Table 4.2 (1, 0, 12, 12, 0, 11) has an adjusted effective size of
#' 3.11, printed 3.17 (the unadjusted value), and the lower end of `E'` in
#' Table 6.6 is 1 / (1 - 0.7015) = 3.350, printed 3.503.
#'
#' @param choice The lineup position each mock witness picked: whole
#'   numbers from 0, 0 for a mock witness who chose no one; one member at
#'   least must be chosen.
#' @param k `NULL` for the largest position picked, or the lineup's size:
#'   one whole number from 2 to 67108864 (2^26).
#' @param susp The suspect's position: one whole number from 1 to `k`.
#' @inheritParams lnroc
#' @return An object of class `lnsize`, a list: `count`, the choices of each
#'   position; `stat`, a named vector (`n`, `k`, `nm` (the mock witnesses
#'   who chose a member), `psusp`, `pbin`, `fsize`, `esize`, `esizeb`,
#'   `etr`, `etrlo`, `etrhi`); and `susp`.
#' @references
#'   Malpass, R. S. (1981). Effective size and defendant bias in eyewitness
#'   identification lineups. *Law and Human Behavior*, 5(4), 299-309.
#'   \doi{10.1007/BF01044945}
#'
#'   Tredoux, C. G. (1996). *Evaluating the Fairness of Identification
#'   Parades with Measures of Facial Similarity*. PhD thesis, University of
#'   Cape Town.
#'
#'   Tredoux, C. G. (1998). Statistical inference on measures of lineup
#'   fairness. *Law and Human Behavior*, 22(2), 217-237.
#'   \doi{10.1023/A:1025746220886}
#'
#'   Tredoux, C. (1999). Statistical considerations when determining
#'   measures of lineup size and lineup bias. *Applied Cognitive
#'   Psychology*, 13(S1), S9-S26.
#'   \doi{10.1002/(SICI)1099-0720(199911)13:1+<S9::AID-ACP634>3.0.CO;2-1}
#'
#'   Agresti, A. and Agresti, B. F. (1978). Statistical analysis of
#'   qualitative variation. *Sociological Methodology*, 9, 204-237.
#'   \doi{10.2307/270810}
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
