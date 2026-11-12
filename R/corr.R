#' Incapacitation: crime rate, share prevented and the marginal year
#'
#' An offender commits crimes at rate `lam` a year while free, and each crime
#' leads to a custodial sentence of `s` years with probability `q`. Over a
#' career the crime rate is `lam / (1 + lam q s)`: the share of offending
#' prevented, `lam q s / (1 + lam q s)`, never reaches 1, longer or surer
#' sentences lower the rate, and one more year of sentence prevents
#' `lam^2 q / ((1 + lam q s) (1 + lam q (s + 1)))` crimes a year, less for
#' each year added. A sentence uniform across offenders prevents a larger
#' share of the offending of those with the higher rate. The model holds
#' the rate fixed over the career and replaces nobody; see [incapr()] for
#' both. See `proof("incap")`.
#'
#' @param lam Offending rates while free (crimes a year): positive numbers,
#'   one per group, optionally named.
#' @param q The chance a crime leads to a custodial sentence: numbers in the
#'   interval from 0 (excluded) to 1, recycled to the length of `lam`.
#' @param s Sentence lengths in years: non-negative numbers, recycled.
#' @param share `NULL`, or one non-negative population share per group,
#'   summing to 1.
#' @return A data frame with one row per group: `group`, `lam`, `q`, `s`,
#'   `rate` (the long-run crime rate), `prev` (the share prevented) and
#'   `marg` (the crimes a year one more year of sentence prevents). With
#'   `share`, a last row `all` holds the weighted free rate in `lam`, the
#'   weighted rate in `rate` and the share of all offending prevented in
#'   `prev`.
#' @references Avi-Itzhak, B. and Shinnar, R. (1973). Quantitative models in
#'   crime control. *Journal of Criminal Justice*, 1(3), 185-217.
#'   \doi{10.1016/0047-2352(73)90061-5}
#' @examples
#' incap(lam = c(low = 2, high = 10), q = 0.1, s = 1, share = c(0.8, 0.2))
#' @export
incap <- function(lam, q, s, share = NULL) {
  .Call(C_clab_r_incap, lam, q, s, share)
}

#' Incapacitation under desistance and replacement
#'
#' Crimes prevented by a sentence of `s` periods served from career age
#' `t0`, when the offending rate falls with career age and a share `r` of
#' the crimes prevented is committed by someone else instead. The count is
#' at most `s` times the rate at the start of the sentence and at least `s`
#' times the rate at its end; the same sentence served one period later
#' prevents no more; and, net of replacement, the estimate from the rate at
#' the start is still an upper bound. The usual constant-rate estimate
#' therefore overstates what a sentence buys under both assumptions. See
#' `proof("incapr")`.
#'
#' @param lam Offending rates by career age, `lam[1]` at age 0: non-negative
#'   numbers that never increase.
#' @param t0 The career age at which the sentence starts: a whole number
#'   from 0, with `t0 + s` at most `length(lam) - 1`.
#' @param s The sentence in periods: a whole number from 1.
#' @param r The replacement share: one number from 0 to 1.
#' @return A named numeric vector: `prev` (crimes prevented), `net` (net of
#'   replacement), `upper` (`s * lam` at age `t0`), `lower` (`s * lam` at age
#'   `t0 + s`) and `later` (crimes prevented by the same sentence from age
#'   `t0 + 1`).
#' @references Shinnar, S. and Shinnar, R. (1975). The effects of the
#'   criminal justice system on the control of crime: a quantitative
#'   approach. *Law and Society Review*, 9(4), 581-611.
#'   \doi{10.2307/3053340}
#' @examples
#' incapr(lam = c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1), t0 = 2, s = 3, r = 0.25)
#' @export
incapr <- function(lam, t0, s, r = 0) {
  .Call(C_clab_r_incapr, lam, t0, s, r)
}

#' Evaluating a release rule from the released alone
#'
#' Failures (a missed court date, a new arrest) are seen only for the
#' defendants a judge released. A proposed rule that releases only people
#' the judge released is scored exactly from their outcomes; one that
#' releases people the judge detained has a failure rate known only up to
#' the worst and best cases for them, an interval whose width is the share
#' of the rule's releases the judge detained. Both ends can occur. The
#' naive score on the rule's observed releases alone always lies inside the
#' interval, but reports one point of it as if it were known. See
#' `proof("sellab")`.
#'
#' @param y Failures: 1 or 0 for each released defendant; values for
#'   detained defendants are ignored and may be `NA`.
#' @param rel Released by the judge: `TRUE` or `FALSE` for each defendant.
#' @param rule Released by the rule: `TRUE` or `FALSE`, same length.
#' @param w `NULL` for equal weights, or one non-negative weight per
#'   defendant.
#' @return A named numeric vector: `obs` (the failure rate among the judge's
#'   releases), `unobs` (the share of the rule's releases the judge
#'   detained), `lower` and `upper` (the bounds on the rule's failure rate),
#'   `width` (`upper - lower`, equal to `unobs`) and `naive` (the failure
#'   rate on the rule's releases that were observed; `NaN` when none were).
#' @references Lakkaraju, H., Kleinberg, J., Leskovec, J., Ludwig, J. and
#'   Mullainathan, S. (2017). The selective labels problem. *Proceedings of
#'   the 23rd ACM SIGKDD International Conference on Knowledge Discovery and
#'   Data Mining*, 275-284. \doi{10.1145/3097983.3098066}
#' @examples
#' set.seed(4)
#' risk <- runif(500)
#' y <- rbinom(500, 1, risk)
#' sellab(y, rel = risk < 0.7, rule = risk < 0.8)
#' @export
sellab <- function(y, rel, rule, w = NULL) {
  .Call(C_clab_r_sellab, y, rel, rule, w)
}

#' Mean time in custody from completed and current stays
#'
#' The mean length of stay of the people released (or of the cases
#' disposed) leaves out everyone still inside. Their elapsed time is a
#' lower bound on their stay, so the mean over everyone admitted is at
#' least `lower`, the mean of the completed stays and the elapsed times
#' together; the released-only mean falls short of it by
#' `m / (n + m)` times the excess of the mean elapsed time over the
#' released mean, and understates the truth for certain when that excess is
#' positive. No upper bound exists without a cap on how long the current
#' stays will last. See `proof("losbnd")`.
#'
#' @param done Completed stays: non-negative numbers.
#' @param open Elapsed time of the stays still open: non-negative numbers.
#' @return A named numeric vector: `mean` (of the completed stays), `age`
#'   (the mean elapsed time of the open ones), `lower` (the lower bound on
#'   the mean over all) and `bias` (`lower - mean`, the least by which the
#'   completed mean understates when it is positive).
#' @references Manski, C. F. (2003). *Partial Identification of Probability
#'   Distributions*. Springer. \doi{10.1007/b97478}
#' @examples
#' losbnd(done = c(30, 45, 60, 90, 120), open = c(100, 150, 200))
#' @export
losbnd <- function(done, open) {
  .Call(C_clab_r_losbnd, done, open)
}

#' Bounds on a sentence's effect without assumptions
#'
#' Each person receives one of two sentences and shows the outcome (say,
#' reconviction) under that one only. With no assumption about who receives
#' which, the chance of the outcome under a sentence is known only up to an
#' interval whose width is the share who received the other sentence, and
#' the effect of one sentence against the other lies in an interval of
#' width exactly 1 that always contains 0: observational data alone cannot
#' even sign it. See `proof("sntbnd")`.
#'
#' @param y Outcomes: 0 or 1 for each person.
#' @param z Sentences: a vector with exactly two distinct values, same
#'   length: character, numbers, logical or a factor (its labels; level order
#'   and whether it is ordered are ignored), not dates.
#' @param treat The value of `z` taken as the treatment, matched exactly;
#'   `NULL` for the second in sorted order. The values are sorted as text
#'   by Unicode code point (`"B"` before `"a"`, `"10"` before `"2"`), the
#'   same on every machine whatever its locale.
#' @param w `NULL` for equal weights, or one non-negative weight per person
#'   with a positive total.
#' @return A data frame with three rows, the comparison sentence, the
#'   treatment and `effect` (treatment minus comparison): `term`, `share`
#'   (the weighted share sentenced; `NA` for the effect), `lower`, `upper`
#'   and `observed` (the observed rate under each sentence and their
#'   difference).
#' @references Manski, C. F. and Nagin, D. S. (1998). Bounding disagreements
#'   about treatment effects: a case study of sentencing and recidivism.
#'   *Sociological Methodology*, 28(1), 99-137.
#'   \doi{10.1111/0081-1750.00043}
#' @examples
#' set.seed(1)
#' z <- sample(c("community", "custody"), 500, TRUE, c(0.6, 0.4))
#' y <- rbinom(500, 1, ifelse(z == "custody", 0.55, 0.35))
#' sntbnd(y, z)
#' @export
sntbnd <- function(y, z, treat = NULL, w = NULL) {
  .Call(C_clab_r_sntbnd, y, z, treat, w)
}

#' Period-two hazards after selective depletion
#'
#' A share `s` of a cohort is at high risk of reconviction, `h` a period,
#' and the rest at low risk, `l`. If a programme changes who survives the
#' first period without reconviction, the period-two hazard among the
#' survivors is a weighted average of `h` and `l` with weight the share of
#' high-risk survivors, and it rises with that share. An arm that keeps
#' more of its high-risk people out of reconviction in period one therefore
#' shows the higher period-two hazard even if it has no effect on anyone in
#' period two: a period-by-period hazard ratio is not a causal contrast; a
#' risk difference at a fixed follow-up is. See `proof("hzdsel")`.
#'
#' @param s The high-risk share at the start, strictly between 0 and 1.
#' @param h The period hazard of the high-risk type: at most 1.
#' @param l The period hazard of the low-risk type: from 0, below `h`.
#' @param surh Period-one survival of the high-risk type in the control and
#'   treated arms: two numbers above 0 and at most 1.
#' @param surl Period-one survival of the low-risk type, the same way.
#' @return A named numeric vector: `share0`, `share1` (the high-risk share
#'   among the survivors of each arm), `haz0`, `haz1` (the period-two hazard
#'   in each arm) and `ratio` (`haz1 / haz0`).
#' @references Hernan, M. A. (2010). The hazards of hazard ratios.
#'   *Epidemiology*, 21(1), 13-15. \doi{10.1097/EDE.0b013e3181c1ea43}
#' @examples
#' hzdsel(s = 0.5, h = 0.5, l = 0.1, surh = c(0.5, 0.8), surl = c(0.9, 0.9))
#' @export
hzdsel <- function(s, h, l, surh, surl) {
  .Call(C_clab_r_hzdsel, s, h, l, surh, surl)
}
