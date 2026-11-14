#' Differences between groups in a risk tool's decisions
#'
#' How the share of cases a tool flags (predicts positive) differs between
#' groups, as a difference and as a ratio against a base group, each with
#' an interval; with the outcomes, how its true and false positive rates
#' differ; and the bias amplification score, the largest parity gap times
#' the Gini coefficient of the groups' shares.
#'
#' `fair()` reports each group's rates and their ratios; `fairdf()` reports
#' the comparisons the fairness literature names. The parity gap `dp` is a
#' group's share flagged minus the base group's (demographic parity holds
#' when every gap is 0), with Newcombe's (1998) hybrid score interval. The
#' ratio `di` is the disparate impact ratio, with the log interval of Katz
#' et al. (1978); under the four-fifths rule a ratio below 0.8 against the
#' group with the highest selection rate is the usual sign of adverse
#' impact, so by default that group is the base. With the outcomes `y`,
#' `dtpr` and `dfpr` are the gaps in true and false positive rates (both 0
#' is equalized odds), `aod` their mean (the average odds difference) with
#' the interval of Zou and Donner (2008) built from the two Newcombe
#' intervals, and `eo` the larger of the two gaps in size. The bias
#' amplification score of Barman and Barman (2026) is large only when a
#' sizeable gap comes with unequal shares across all the groups; its
#' bootstrap interval, and the Gini coefficient's, are percentile intervals
#' from redrawing each group's decisions with replacement (Efron and
#' Tibshirani 1993), the base group held.
#'
#' These numbers describe what the tool did to each group. They do not say
#' why: a gap in shares flagged can come from a gap in the outcomes the
#' tool was built to predict, from how those outcomes were recorded, or
#' from the tool; and when the groups' base rates differ, equal error rates
#' and equal predictive values cannot all hold at once (Chouldechova 2017,
#' see [impfpr()]). A ratio's interval is `NaN` when either share is 0, and
#' a ratio against a base group with no case flagged is `Inf` (or `NaN` for
#' 0 over 0).
#'
#' @param pred The tool's decisions: 1 (or `TRUE`) for a case it flags, 0
#'   (or `FALSE`) otherwise, one per case.
#'   A `clpre` list from [clprep()] is taken too: its elements fill the
#'   arguments by name.
#' @param y `NULL`, or the outcomes: 0 or 1 (or `FALSE`/`TRUE`), one per
#'   case.
#' @param base The base group: `NULL` for the group with the largest share
#'   flagged (the first in order on ties), its label, or its position from
#'   1.
#' @param nboot The bootstrap replicates for the intervals of the Gini
#'   coefficient and the score: a whole number from 0 (none).
#' @param seed `NULL` to draw the seed of the bootstrap (or the permutations)
#'   from R's generator, so that [set.seed()] fixes it, or one whole number
#'   from 0.
#' @inheritParams fair
#' @return A list: `tab`, a data frame with one row per group (in the order
#'   of `grp`): `group`, `n`, `npos` (cases flagged), `ppos` (their share)
#'   with `lo` and `hi`, `dp` with `dplo` and `dphi`, `di` with `dilo` and
#'   `dihi`, and with `y` also `tpr`, `fpr`, `dtpr` (`dtprlo`, `dtprhi`),
#'   `dfpr` (`dfprlo`, `dfprhi`), `aod` (`aodlo`, `aodhi`) and `eo`; the base
#'   group's own intervals are `NaN`. `sum`, a named vector: `di` (the
#'   smallest finite ratio of another group), `dp` (the gap of largest size,
#'   signed), `gini`, `bas`, `ginilo`, `ginihi`, `baslo`, `bashi` (`NaN` when
#'   `nboot` is 0) and with `y` also `aod` and `eo` (those of largest size,
#'   signed). `base`, the base group's label.
#' @references Newcombe, R. G. (1998). Interval estimation for the
#'   difference between independent proportions: comparison of eleven
#'   methods. *Statistics in Medicine*, 17(8), 873-890.
#'   \doi{10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I}
#'
#'   Katz, D., Baptista, J., Azen, S. P. and Pike, M. C. (1978). Obtaining
#'   confidence intervals for the risk ratio in cohort studies.
#'   *Biometrics*, 34(3), 469-474. \doi{10.2307/2530610}
#'
#'   Zou, G. Y. and Donner, A. (2008). Construction of confidence limits
#'   about effect measures: a general approach. *Statistics in Medicine*,
#'   27(10), 1693-1702. \doi{10.1002/sim.3095}
#'
#'   Barman, P. K. and Barman, P. K. (2026). Spatial disparities in
#'   simulated patrol detection: a reproducible generative-model audit
#'   across Baltimore and Chicago. \doi{10.48550/arXiv.2603.18987}
#'
#'   Efron, B. and Tibshirani, R. J. (1993). *An Introduction to the
#'   Bootstrap*. Chapman and Hall. \doi{10.1007/978-1-4899-4541-9}
#' @examples
#' set.seed(1)
#' g <- rep(c("a", "b"), each = 200)
#' y <- rbinom(400, 1, 0.3)
#' p <- rbinom(400, 1, ifelse(g == "a", 0.4, 0.25))
#' r <- fairdf(p, g, y, nboot = 199, seed = 1)
#' r$tab[, c("group", "ppos", "dp", "dplo", "dphi", "di", "dilo", "dihi")]
#' r$sum
#' @export
fairdf <- function(pred, grp = NULL, y = NULL, base = NULL, conf = 0.95,
                   nboot = 0, seed = NULL) {
  if (inherits(pred, "clpre")) {
    return(.prep(fairdf, pred, match.call(), parent.frame()))
  }
  .Call(C_clab_r_fairdf, pred, grp, y, base, conf, nboot, seed)
}

#' Disparity metrics over time and place
#'
#' The disparate impact ratio, the parity gap, the Gini coefficient of the
#' groups' shares flagged and the bias amplification score of [fairdf()],
#' for every city and period, and how far each city's ratio moves over the
#' periods (after Barman and Barman 2026).
#'
#' A cell (one city in one period) is audited when it holds at least two
#' groups and the base group; otherwise it is reported and skipped. The
#' base group is the same in every cell: by default the group with the
#' largest share flagged over all the data. For each city, the means of the
#' four values over its audited cells (finite values only), the smallest
#' and largest ratio, their range, and how many periods have a ratio above
#' 1. A ratio that moves a lot from one period to the next says that one
#' audit at one time describes that time only; a spread between cities
#' says the same tool can treat groups differently in different places.
#' Neither says why. For intervals on one cell, call [fairdf()] on it.
#'
#' @param period The period of each case: labels (text, numbers or a
#'   factor), sorted as [fair()] sorts groups; ISO dates as text sort in
#'   time order.
#' @param city `NULL` (one place, labelled `"all"`), or the city of each
#'   case: labels.
#' @inheritParams fairdf
#' @return A list: `cells`, a data frame with one row per city and period
#'   that occurs: `city`, `period`, `n`, `ok` (`TRUE` when audited), `di`,
#'   `dp`, `gini` and `bas` (`NaN` when skipped); `cities`, one row per
#'   city: `city`, `ncell` (cells audited), `di`, `dp`, `gini`, `bas` (their
#'   means), `dimin`, `dimax`, `dirng` and `ngt1`; `sum`, `rng` (the largest
#'   range) and `spread` (largest minus smallest mean ratio of the cities;
#'   `NaN` with fewer than two); `base`.
#' @references
#'   Barman, P. K. and Barman, P. K. (2026). Spatial disparities in
#'   simulated patrol detection: a reproducible generative-model audit
#'   across Baltimore and Chicago. \doi{10.48550/arXiv.2603.18987}
#' @examples
#' set.seed(2)
#' n <- 600
#' g <- sample(c("a", "b"), n, TRUE)
#' per <- sample(c("2024-1", "2024-2", "2024-3"), n, TRUE)
#' cty <- sample(c("north", "south"), n, TRUE)
#' p <- rbinom(n, 1, ifelse(g == "a", 0.4, ifelse(cty == "north", 0.3, 0.2)))
#' tmpaud(p, g, per, cty)$cities
#' @export
tmpaud <- function(pred, grp = NULL, period = NULL, city = NULL,
                   base = NULL) {
  if (inherits(pred, "clpre")) {
    return(.prep(tmpaud, pred, match.call(), parent.frame()))
  }
  .Call(C_clab_r_tmpaud, pred, grp, period, city, base)
}

#' Individual fairness: a Lipschitz audit
#'
#' Dwork et al. (2012) ask that similar individuals be treated similarly:
#' that the gap between two people's scores be at most `lip` times the
#' distance between them, for a distance that says who is alike for the
#' task. `lipaud()` checks the given pairs against that rule and reports
#' the smallest constant the scores satisfy on them.
#'
#' A pair at distance 0 is two people the distance cannot tell apart: any
#' gap between their scores breaks the rule, and no constant repairs it
#' (its ratio is `Inf`). The audit is only as good as the distance: it says
#' whether the scores treat alike those the distance calls alike, not
#' whether the distance is the right one.
#'
#' @param h The scores: finite numbers, one per person.
#' @param d The distances: with `pairs` `NULL`, a [stats::dist()] object or
#'   a symmetric matrix between all the people (every pair is audited);
#'   with `pairs`, `NULL` (all 1) or one number from 0 per pair.
#' @param pairs `NULL`, or the pairs to audit: a matrix of two columns of
#'   person numbers from 1.
#' @param lip The constant: one number from 0.
#' @return A list: `lreq` (the smallest constant the scores satisfy on these
#'   pairs), `nviol` and `rate` (the pairs breaking the rule at `lip`, and
#'   their share), `ok` (`TRUE` when none does), `maxgap` (the largest
#'   gap), `worst` (the first pair with the largest ratio), and `pairs`, a
#'   data frame with one row per pair: `i`, `j`, `d`, `gap`, `ratio` and
#'   `viol`.
#' @references Dwork, C., Hardt, M., Pitassi, T., Reingold, O. and Zemel,
#'   R. (2012). Fairness through awareness. *Proceedings of the 3rd
#'   Innovations in Theoretical Computer Science Conference*, 214-226.
#'   \doi{10.1145/2090236.2090255}
#' @examples
#' x <- c(1, 2, 2.5, 6)
#' h <- c(0.10, 0.20, 0.60, 0.65)
#' # scores on a scale of 0 to 1 against distances in years
#' lipaud(h, dist(x), lip = 0.2)
#' @export
lipaud <- function(h, d = NULL, pairs = NULL, lip = 1) {
  .Call(C_clab_r_lipaud, h, d, pairs, lip)
}
