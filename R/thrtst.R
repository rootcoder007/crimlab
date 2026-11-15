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
#'   With `fit = "ml"` (the default) the model is fitted by maximum
#'   likelihood to the binomial counts of searches among stops and hits
#'   among searches. A lower threshold for a group is
#'   evidence of discrimination against it. `thr` is each group's threshold
#'   averaged over the areas, weighted by the areas' stops, with a
#'   delta-method standard error; the likelihood-ratio test compares the
#'   fit with one threshold for every group. The model is fitted only when
#'   it has more rates (two per group and area with stops) than parameters
#'   (`3 (groups + areas - 1)`): with `area` `NULL` it is not.
#' * The tests of Anwar and Fang (2006), when `off` gives the officer's
#'   group. If officers of no group are prejudiced relative to the others,
#'   the ranking of the officer groups by search rate is the same for
#'   people of every group, and likewise the ranking by hit rate
#'   (Proposition 3; Section 2.1.2). `mono` is their test of monolithic
#'   behaviour (Proposition 2, equations 3 and 4; Section 3.2.3): for each
#'   group of the people stopped and each measure, Pearson's chi-square of
#'   equal rates across the officer groups, on (officer groups - 1) df.
#'   With `rank = "pair"` (the default) each pair of officer groups is
#'   taken alone: for each measure the z statistics of the difference
#'   across the suspects' groups are given, and `rej` is `TRUE` when one
#'   is significantly positive and another significantly negative at
#'   `level` (two-sided), a ranking that changes with the suspect's group.
#'   With `rank = "joint"` the ranking of all officer groups is tested at
#'   once, as the paper states it: each one-sided z of their equation 7
#'   (the difference of two rates over the root of the sum of the sample
#'   variances of the 0/1 outcomes over the counts) beyond the normal
#'   quantile at `level` orders two officer groups among the people of one
#'   group; the null of no prejudice holds when one ranking agrees with
#'   every such order in every group, and is rejected (`rej`) when the
#'   orders form a cycle (a above b in one group, b above c in another, c
#'   above a in a third, or a reversed pair). The common `rank` (1 the
#'   highest rate, ties broken by the pooled rate) is given when it is not
#'   rejected. As the paper notes, the test has power below one: a group
#'   of officers prejudiced against every group alike is not detected. The
#'   paper's p-values come from the normal and chi-square distributions,
#'   as here; its resampling (Section 2.2) that balances the officers'
#'   assignment across stations needs data on single officers and is not
#'   done.
#'
#' With `fit = "bayes"` the threshold test is fitted as the papers fit it,
#' by Bayesian inference. The signal distributions are those above, with
#' the area with the most stops as the reference (its terms 0), and each
#' group and area has its own threshold. The priors are those of Simoiu,
#' Corbett-Davies and Goel (2017, Section 2): normal priors with standard
#' deviation 2 on the groups' terms `phi_r` and `delta_r`; the areas'
#' terms drawn from `N(mu, sigma)`, with `mu ~ N(0, 2)` and `sigma` half
#' normal with scale 2, for each of `phi` and `delta`; and each threshold
#' on the logit scale drawn from `N(mu_t, sigma_t)` of its group, again
#' with `mu_t ~ N(0, 2)` and half-normal `sigma_t`. Pierson,
#' Corbett-Davies and Goel (2017) keep this model and put the
#' discriminant distributions in place of the beta distributions of
#' Simoiu et al. The posterior is sampled by Hamiltonian Monte Carlo, as
#' both papers do with Stan, here with the No-U-Turn sampler of Hoffman
#' and Gelman (2014) and a diagonal mass matrix: `chains` chains of `iter`
#' steps, the first half of each tuning the step size (to an acceptance
#' of 0.8) and the scale of each parameter and then dropped, at most 2^10
#' leapfrog steps a draw; `seed` fixes the draws. Each chain starts at the
#' maximum-likelihood fit with every parameter moved at random by up to
#' 0.25, because random starts can leave a chain in a far worse mode.
#' `thr` then holds each group's threshold averaged over its areas by
#' their stops (posterior mean `thr`, standard deviation `se` and the
#' interval `lo`, `hi` of probability `level`) and its difference from
#' the reference group's (`diff`, `dse`, `dlo`, `dhi`); `post` gives the
#' same summaries with the
#' effective sample size and the split R-hat of Gelman et al. (2013) for
#' every parameter. The function refuses to report chains that have not
#' mixed: an R-hat of 1.05 or more for any parameter (the rule of Simoiu
#' et al., Section 3.1), or fewer than 100 effective draws of a threshold
#' or a difference; more iterations then help. The data identify the
#' model only when there are more cells with stops than `2 (groups +
#' areas - 1)`, three groups and five areas at least (Simoiu et al.,
#' Section 2); with fewer the function refuses. The authors' published
#' code for the faster model sets other priors (a threshold scale shared
#' by the groups, tighter priors on the areas); the priors here are those
#' the paper states.
#'
#' The two fits answer under different models and should be reported as
#' such, not as checks of one another. The Bayesian model gives each
#' group and area its own threshold under a hierarchical prior; on data
#' made with one threshold per group (and area offsets) it can report
#' thresholds far from the maximum-likelihood ones and from the values
#' that made the data. For three groups in ten areas with 20,000 stops a
#' cell its posterior means were about 0.13, 0.07 and 0.14, against 0.27,
#' 0.20 and 0.31 from maximum likelihood and from the generating model.
#' Stan's sampler gives the same posterior on the same model, so this is
#' a property of the model, not of the sampler.
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
#' @param rank `"pair"` or `"joint"`: the rank-order test by pairs of officer
#'   groups or for all of them at once.
#' @param fit `"ml"` or `"bayes"`: the threshold model by maximum likelihood
#'   or by Bayesian inference.
#' @param chains The chains of the Bayesian fit: one whole number from 2.
#' @param iter The steps of each chain, half of them warm-up: one whole
#'   number from 20.
#' @param seed The seed of the Bayesian fit: `NULL` to draw it from R's
#'   generator, or one whole number from 0.
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
#'   (`model`, `conv`, `iter`; with `fit = "bayes"` one row with `chains`,
#'   the divergent transitions `div`, the mean acceptance `acc` and the
#'   step size `eps` as well); `cells`, the fitted `phi`, `delta`, `thr`,
#'   `srate` and `hrate` of each group and area with stops; `rank`, `NULL`
#'   without `off`, otherwise with `rank = "pair"` one row per pair of
#'   officer groups and measure (`off1`, `off2`, `meas`, `zmin`, `zmax`,
#'   `rej`) and with `rank = "joint"` one row per measure and officer group
#'   (`meas`, `off`, `rank`, `rej`); `ref`; and `mono`, `NULL` without
#'   `off`, otherwise one row per group of the people stopped and measure
#'   (`grp`, `meas`, `chisq`, `df`, `p`); `post`, `NULL` unless `fit =
#'   "bayes"`, otherwise one row per quantity (`param`, `mean`, `sd`, `lo`,
#'   `hi`, `ess`, `rhat`). With `fit = "bayes"` `thr` has the columns
#'   `group`, `thr`, `se`, `lo`, `hi`, `diff`, `dse`, `dlo` and `dhi`, and
#'   `lrt` is `NULL`.
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
#'
#'   Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A.
#'   and Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd ed. Chapman and
#'   Hall/CRC. \doi{10.1201/b16018}
#'
#'   Hoffman, M. D. and Gelman, A. (2014). The No-U-Turn sampler:
#'   adaptively setting path lengths in Hamiltonian Monte Carlo. *Journal
#'   of Machine Learning Research*, 15, 1593-1623.
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
#' # the joint test of three officer groups
#' thrtst(rep(600, 9), c(72, 48, 60, 60, 72, 48, 48, 60, 72),
#'        rep(c(36, 24, 30), 3), grp = rep(1:3, each = 3),
#'        off = rep(c("a", "b", "c"), 3), rank = "joint")$rank
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
  maxit = 2000,
  rank = "pair",
  fit = "ml",
  chains = 4,
  iter = 1000,
  seed = NULL
) {
  .Call(
    C_clab_r_thrtst, stops, srch, hits, grp, area, off, ref, level, maxit,
    rank, fit, chains, iter, seed
  )
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
