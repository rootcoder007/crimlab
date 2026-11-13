#' Calibration audit of a predictive-policing ranking
#'
#' Does the order in which an algorithm ranks areas by predicted risk match
#' the order of what was later recorded there, and does the disagreement
#' follow the areas' demographic make-up? `calaud()` ranks the areas both
#' ways, measures the agreement by Spearman's correlation, and gives each
#' group the mean gap between its areas' two ranks, with a permutation
#' test. It generalises the district-level analysis of the SciencesPo
#' predictive-policing study of Chicago (Lacherade et al. 2021) to any
#' place: the five columns of the audit are `area`, `risk`, `outcome`,
#' `group` and `population`.
#'
#' Records are rolled up to areas first: an area's risk is the mean over
#' its records, its rate the mean outcome (or, with `population`, the
#' summed outcome per head), and its group the most frequent label among
#' its records (the first in order on ties). One record per area is the
#' per-area audit itself. Rank 1 is the highest; tied values share the mean
#' of their ranks. An area's gap is its outcome rank minus its risk rank:
#' positive when the algorithm ranks it as riskier than its recorded
#' outcome warrants, negative when less risky.
#'
#' Spearman's rho is Pearson's correlation of the two rankings (Spearman
#' 1904). `pex` is its two-sided p-value by algorithm AS 89 (Best and
#' Roberts 1975): exact for nine areas or fewer, from the Edgeworth series
#' up to 1289, and `NA` when either ranking has ties; `pt` is the
#' large-sample one from the t distribution, as [stats::cor.test()] gives
#' them. A group's mean gap is tested against shuffling the groups across
#' the areas, `nperm` times; `between` tests all groups at once by the
#' between-group sum of squares of the gaps (the gaps sum to 0).
#'
#' The audit compares the algorithm with the records, not with crime: a
#' recorded outcome follows where the police looked as well as what
#' happened, so a group's areas can look over-predicted because their
#' incidents were recorded less, or under-predicted because they were
#' recorded more (see [fairsim()]). A gap says the two orders disagree for
#' that group's areas; it does not say which order is right or why.
#'
#' @param area The area of each record: labels; or a data frame with
#'   columns `area`, `risk`, `outcome`, `group` and optionally
#'   `population`, the other arguments then left out.
#'   A `clpre` list from [clprep()] is taken too: its elements fill the
#'   arguments by name.
#' @param risk The risk the algorithm predicts for each record: finite
#'   numbers.
#' @param outcome The recorded outcome of each record (a 0 or 1, or a
#'   count): finite numbers.
#' @param group The group of each record (or area): labels.
#' @param population `NULL`, or the population of each area: positive, one
#'   per record (the same within an area) or one per area named by area.
#' @param nperm The number of permutations: a whole number from 0 (none).
#' @inheritParams fairdf
#' @return A list: `areas`, a data frame with one row per area (sorted by
#'   label): `area`, `n` (records), `risk`, `rate`, `group`, `riskrk`,
#'   `raterk` and `gap`; `groups`, one row per group: `group`, `n` (areas),
#'   `gap` (mean) and `p` (`NaN` without areas or permutations); `spear`, a
#'   named vector: `rho`, `s` (Spearman's statistic, the sum of squared
#'   rank differences), `pex`, `t` and `pt`; `ties` (`TRUE` when either
#'   ranking has ties); `between`, `ss` and `p`; and `nperm`.
#' @references Spearman, C. (1904). The proof and measurement of
#'   association between two things. *American Journal of Psychology*,
#'   15(1), 72-101. \doi{10.2307/1412159}
#'
#'   Best, D. J. and Roberts, D. E. (1975). Algorithm AS 89: the upper tail
#'   probabilities of Spearman's rho. *Journal of the Royal Statistical
#'   Society, Series C*, 24(3), 377-379. \doi{10.2307/2347111}
#'
#'   Lacherade, P., Szabo, B., Krikava, J. and Aeby, M. (2021).
#'   Predictive policing in the City of Chicago. Group project, Decoding
#'   Biases in AI, Sciences Po.
#'   <https://github.com/marquito95/Predictive-policing-Chicago>
#' @examples
#' set.seed(3)
#' d <- fairsim(4000, nareas = 24, bias = 0.5, seed = 1)
#' # the recorded outcome is `det`: the incidents that were detected
#' r <- calaud(d$area, d$risk, d$det, d$group, nperm = 499, seed = 1)
#' r$groups
#' r$spear
#' @export
calaud <- function(area, risk = NULL, outcome = NULL, group = NULL,
                   population = NULL, nperm = 999, seed = NULL) {
  if (inherits(area, "clpre")) {
    return(.prep(calaud, area, match.call(), parent.frame()))
  }
  .Call(C_clab_r_calaud, area, risk, outcome, group, population, nperm, seed)
}

#' A risk score's spread across groups
#'
#' Each group's number of cases, mean, median and standard deviation of a
#' score, each mean's gap to a reference group with its interval, and the
#' one-way analysis of variance of the score on group, with Welch's (1951)
#' test, which does not assume equal variances.
#'
#' The analysis of variance is the standard one ([stats::anova()] of
#' [stats::lm()]): a group with one case adds to the between-group sum of
#' squares and nothing to the within-group one. The gaps' intervals use the
#' pooled within-group variance on its degrees of freedom. Welch's test
#' needs two or more cases and some spread in every group, and is `NaN`
#' otherwise.
#'
#' A gap in mean scores is not by itself a sign of bias: a score built to
#' predict an outcome differs between groups whose rates of that outcome
#' differ. Read it with [calaud()], which compares the scores with what was
#' recorded, or [fair()], which compares them with outcomes case by case.
#'
#' @param score The scores: finite numbers, one per case.
#'   A `clpre` list from [clprep()] is taken too: its elements fill the
#'   arguments by name.
#' @param group The group of each case: labels.
#' @param ref The reference group: `NULL` for the group with the lowest mean
#'   (the first on ties), its label, or its position from 1.
#' @inheritParams fair
#' @return A list: `tab`, a data frame with one row per group: `group`, `n`,
#'   `mean`, `median`, `sd` (`NaN` for one case), `gap`, `se`, `lo` and `hi`
#'   (`NaN` for the reference); `anova` and `welch`, named vectors `f`,
#'   `df1`, `df2` and `p`; `spread` (the largest mean minus the smallest);
#'   `ref`.
#' @references Welch, B. L. (1951). On the comparison of several mean
#'   values: an alternative approach. *Biometrika*, 38(3-4), 330-336.
#'   \doi{10.1093/biomet/38.3-4.330}
#' @examples
#' set.seed(7)
#' g <- rep(c("x", "y", "z"), c(40, 40, 20))
#' s <- rnorm(100, ifelse(g == "x", 0, ifelse(g == "y", 0.5, 1)))
#' scrdsp(s, g)
#' @export
scrdsp <- function(score, group = NULL, ref = NULL, conf = 0.95) {
  if (inherits(score, "clpre")) {
    return(.prep(scrdsp, score, match.call(), parent.frame()))
  }
  .Call(C_clab_r_scrdsp, score, group, ref, conf)
}
