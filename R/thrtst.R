#' Tests of discrimination in searches: hit rates, thresholds, rank order
#'
#' Asks whether the police search people of one group on less evidence
#' than people of another. Three tests read the same counts of stops,
#' searches and successful searches ("hits", contraband found):
#'
#' * The outcome (hit-rate) test of Knowles, Persico and Todd (2001): if
#'   officers search whenever the chance of finding something is high
#'   enough, a group searched on less evidence shows a lower hit rate. The
#'   hit rate of each group is compared with the reference group's
#'   (difference, Wald standard error, z and two-sided p) and all groups
#'   together by Pearson's chi-square of hits and misses among the
#'   searched. A hit rate can differ between groups that face the same
#'   threshold when their risks are spread differently (the problem of
#'   infra-marginality), which the threshold test addresses.
#' * The threshold test of Simoiu, Corbett-Davies and Goel (2017), in the
#'   fast form of Pierson, Corbett-Davies and Goel (2017). Each stopped
#'   person carries a signal, normal with mean `delta` for those who carry
#'   contraband (a share `phi` of the stopped) and mean 0 for the rest; the
#'   officer's estimate of the chance of a hit is the posterior given the
#'   signal (a "discriminant distribution"), and a search follows when that
#'   chance exceeds the threshold `t`. In each group and area `phi`,
#'   `delta` and `t` are the logistic (`phi`, `t`) or exponential (`delta`)
#'   of a group term plus an area term, so a group's threshold may differ
#'   from another's by the same amount on the logit scale in every area.
#'   The model is fitted by maximum likelihood to the binomial counts of
#'   searches among stops and hits among searches (where the published
#'   method uses Bayesian inference). A lower threshold for a group is
#'   evidence of discrimination against it. `thr` is each group's threshold
#'   averaged over the areas, weighted by the areas' stops, with a
#'   delta-method standard error; the likelihood-ratio test compares the
#'   fit with one threshold for every group. The model is fitted only when
#'   it has more rates (two per group and area with stops) than parameters
#'   (`3 (groups + areas - 1)`): with `area` `NULL` it is not.
#' * The rank-order test of Anwar and Fang (2006), when `off` gives the
#'   officer's group: if officers of no group are prejudiced relative to
#'   the others, officers of one group search people of every group more
#'   often than officers of another (or less often, in every group), and
#'   likewise for hit rates. For each pair of officer groups and each
#'   measure the z statistics of the difference across the suspects'
#'   groups are given; `rej` is `TRUE` when one is significantly positive
#'   and another significantly negative at `level`, a ranking that changes
#'   with the suspect's group. This is a pairwise simplification of their
#'   test, not its joint form.
#'
#' Rows are cells of a table (or single stops, with counts of 0 and 1):
#' rows of the same group and area (and officer group) add up.
#'
#' @param stops The stops of each row: numbers from 0.
#' @param srch The searches among them: numbers from 0, at most `stops`.
#' @param hits The successful searches: numbers from 0, at most `srch`.
#' @param grp The group of the people stopped: labels, two groups at least,
#'   each with searches.
#' @param area `NULL` or the area (department, precinct) of each row:
#'   labels.
#' @param off `NULL` or the group of the officers of each row: labels, two
#'   groups at least.
#' @param maxit The quasi-Newton steps of each of the three starts of each
#'   fit: one whole number from 1.
#' @inheritParams cntreg
#' @inheritParams kmrec
#' @return An object of class `thrtst`, a list: `hit`, one row per group
#'   (`group`, `stops`, `srch`, `hits`, `srate` and `hrate`, the search and
#'   hit rates, and `diff`, `se`, `z`, `p` of the hit rate against the
#'   reference); `kpt`, the chi-square of equal hit rates (`chisq`, `df`,
#'   `p`); `thr`, one row per group (`group`, `thr`, `se`, `lo`, `hi`, and
#'   `diff`, `dse`, `z`, `p` against the reference), `NULL` when the model
#'   is not identified; `lrt`, the likelihood-ratio test of equal
#'   thresholds (`stat`, `df`, `p`, the log-likelihoods `llfull` and
#'   `llsame`, and the number of `cells`); `fit`, one row per model
#'   (`model`, `conv`, `iter`); `cells`, the fitted `phi`, `delta`, `thr`,
#'   `srate` and `hrate` of each group and area with stops; `rank`, `NULL`
#'   without `off`, otherwise one row per pair of officer groups and
#'   measure (`off1`, `off2`, `meas`, `zmin`, `zmax`, `rej`); and `ref`.
#' @references
#'   Knowles, J., Persico, N. and Todd, P. (2001). Racial bias in motor
#'   vehicle searches: theory and evidence. *Journal of Political
#'   Economy*, 109(1), 203-229. \doi{10.1086/318603}
#'
#'   Anwar, S. and Fang, H. (2006). An alternative test of racial prejudice
#'   in motor vehicle searches: theory and evidence. *American Economic
#'   Review*, 96(1), 127-151. \doi{10.1257/000282806776157579}
#'
#'   Simoiu, C., Corbett-Davies, S. and Goel, S. (2017). The problem of
#'   infra-marginality in outcome tests for discrimination. *The Annals of
#'   Applied Statistics*, 11(3), 1193-1216. \doi{10.1214/17-AOAS1058}
#'
#'   Pierson, E., Corbett-Davies, S. and Goel, S. (2017). Fast threshold
#'   tests for detecting discrimination. arXiv:1702.08536.
#'   \doi{10.48550/arXiv.1702.08536}
#' @examples
#' # stops, searches and hits of two groups in four areas
#' d <- data.frame(
#'   grp = rep(c("a", "b"), 4), area = rep(1:4, each = 2),
#'   stops = c(1200, 800, 900, 1100, 1500, 700, 1000, 1000),
#'   srch = c(60, 70, 50, 95, 80, 66, 52, 90),
#'   hits = c(24, 18, 19, 25, 30, 17, 21, 24)
#' )
#' thrtst(d$stops, d$srch, d$hits, d$grp, area = d$area)
#' # the rank-order test needs the officers' groups
#' thrtst(c(100, 100, 100, 100), c(10, 20, 30, 10), c(5, 5, 5, 5),
#'        grp = c(1, 1, 2, 2), off = c("x", "y", "x", "y"))$rank
#' @export
thrtst <- function(
  stops,
  srch,
  hits,
  grp,
  area = NULL,
  off = NULL,
  ref = NULL,
  level = 0.95,
  maxit = 2000
) {
  .Call(C_clab_r_thrtst, stops, srch, hits, grp, area, off, ref, level, maxit)
}

#' @export
print.thrtst <- function(x, ...) {
  k <- x[["kpt"]]
  cat("thrtst: reference ", x[["ref"]], "\n", sep = "")
  cat(
    "hit rates (chi-square ",
    format(signif(k[["chisq"]], 4)),
    " on ",
    k[["df"]],
    " df, p ",
    format(signif(k[["p"]], 3)),
    "):\n",
    sep = ""
  )
  print(x[["hit"]], row.names = FALSE)
  if (is.null(x[["thr"]])) {
    cat("thresholds: not identified (give `area` with enough areas)\n")
  } else {
    l <- x[["lrt"]]
    cat(
      "thresholds (likelihood ratio ",
      format(signif(l[["stat"]], 4)),
      " on ",
      l[["df"]],
      " df, p ",
      format(signif(l[["p"]], 3)),
      "):\n",
      sep = ""
    )
    print(x[["thr"]], row.names = FALSE)
  }
  if (!is.null(x[["rank"]])) {
    cat("rank order of officer groups:\n")
    print(x[["rank"]], row.names = FALSE)
  }
  invisible(x)
}
