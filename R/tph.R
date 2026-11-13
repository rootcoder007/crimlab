#' Log-ratio transforms of compositions
#'
#' Element shares from portable X-ray fluorescence of bone or soil, or any
#' parts of a whole, live on the simplex: they sum to one, and ordinary
#' statistics on the raw shares mislead. The centred log-ratio (`"clr"`)
#' takes each part's log over the geometric mean of its row; the isometric
#' log-ratio (`"ilr"`, Egozcue et al. 2003) gives `d - 1` coordinates in an
#' orthonormal basis (the default pivot basis: coordinate `i` contrasts the
#' first `i` parts with part `i + 1`), on which the usual methods apply.
#' Each row is closed to sum 1 after `ps` is added to every part, so that a
#' part below the detection limit (0) has a log.
#'
#' @param x The compositions: a numeric matrix or data frame of numbers
#'   from 0, one row per sample and one column per part, at least two
#'   parts.
#' @param type `"clr"` or `"ilr"`.
#' @param ps The pseudocount added to every part: one number from 0.
#' @return A matrix with one row per sample: `d` columns named as the parts
#'   (`"clr"`), or `d - 1` columns `ilr1`, `ilr2`, ... (`"ilr"`).
#' @references
#'   Aitchison, J. (1986). *The Statistical Analysis of Compositional
#'   Data*. Chapman and Hall.
#'
#'   Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and
#'   Barcelo-Vidal, C. (2003). Isometric logratio transformations for
#'   compositional data analysis. *Mathematical Geology*, 35(3), 279-300.
#' @examples
#' x <- rbind(c(Ca = 60, P = 25, Fe = 10, Zn = 5),
#'            c(Ca = 80, P = 12, Fe = 6, Zn = 2))
#' clrilr(x)
#' clrilr(x, "ilr")
#' @export
clrilr <- function(x, type = "clr", ps = 1e-6) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_clrilr, x, type, ps)
}

#' Likelihood ratio of measured evidence
#'
#' Weighs measurements (a tissue density, a pXRF calcium level, a marker of
#' antemortem ingestion) under two competing hypotheses, each a normal
#' model of what the measurement would be: the likelihood ratio is the
#' probability of the evidence under H1 over that under H2. Its verbal band
#' follows the usual forensic scale on `|log10 LR|`: 0 none, up to 1 weak,
#' 2 moderate, 3 moderately strong, 4 strong, 6 very strong, above
#' extremely strong support. The ratio measures support, not the
#' probability that a hypothesis is true.
#'
#' @param x The measurements: finite numbers.
#' @param m1,s1 The mean and standard deviation of a measurement under H1:
#'   one number each or one per measurement (`s1` positive).
#' @param m2,s2 The same under H2.
#' @return A list: `ll1`, `ll2` the log-likelihoods, `llr` their
#'   difference (the natural log of the ratio), `lr` the ratio, `l10` its
#'   log10, and `verbal`, the band and the hypothesis it favours.
#' @references
#'   ENFSI (2015). *ENFSI Guideline for Evaluative Reporting in Forensic
#'   Science*. European Network of Forensic Science Institutes.
#' @examples
#' preslr(c(1.9, 2.3, 2.1), m1 = 2, s1 = 0.3, m2 = 1, s2 = 0.4)
#' @export
preslr <- function(x, m1, s1, m2, s2) {
  .Call(C_clab_r_preslr, x, m1, s1, m2, s2)
}

#' Decay chain of a body
#'
#' A discrete-time Markov chain of decomposition: from each transient stage
#' (fresh, bloat, active and advanced decay, or the user's own stages) a
#' body advances to the next stage with probability `decay * (1 - pres)`
#' (after the last, to skeletal), mummifies with probability `mum * pres`,
#' and otherwise stays.
#' `pres` is the strength of what preserves the body (aridity, quicklime,
#' a sealed casket), 0 for none. From the fundamental matrix come the
#' probability of ending skeletal or mummified from each stage and the
#' expected steps to get there. With no preservation nothing mummifies, so
#' the column `mummified` of `absorb` is also the change in the chance of
#' mummifying that the preservation makes. With `steps` above 0 one path
#' is drawn from `start`.
#'
#' @param pres The preservation: one number from 0 to 1.
#' @param decay The chance to advance a stage without preservation: one
#'   number above 0, at most 1.
#' @param mum The chance to mummify at full preservation: one number from 0
#'   to 1.
#' @param states The names of the transient stages, in order.
#' @param start The stage a drawn path starts from: a name of `states` or
#'   its number.
#' @param steps The longest path drawn (it stops when absorbed): one whole
#'   number from 0, 0 for none.
#' @inheritParams agemix
#' @return A list: `P`, the transition matrix (stages, then `skeletal` and
#'   `mummified`); `N`, the fundamental matrix of the transient stages;
#'   `absorb`, the probability of ending `skeletal` or `mummified` from each
#'   stage; `time`, the expected steps to absorption from each stage;
#'   `start`; `path`, `NULL` or the states of the drawn path.
#' @references
#'   Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*. Van
#'   Nostrand.
#' @examples
#' tpdmc(pres = 0.4)$absorb
#' tpdmc(pres = 0.4, steps = 50, seed = 1)$path
#' @export
tpdmc <- function(pres = 0, decay = 0.5, mum = 0.5,
                  states = c("fresh", "bloat", "active", "advanced"),
                  start = 1, steps = 0, seed = NULL) {
  .Call(C_clab_r_tpdmc, pres, decay, mum, states, start, steps, seed)
}

#' Bayesian preservation model
#'
#' A Bayesian linear model of a preservation score (or any outcome) on
#' burial conditions and measurements, with a normal prior on each
#' coefficient and the noise variance fixed at the least-squares residual
#' variance (an empirical Bayes plug-in), so that the posterior is normal in
#' closed form. An informative prior on a treatment (say, that quicklime
#' preserves) is set through `pm` and `psd`, and the data revise it. With
#' groups (cemeteries, sites) the mean residual of each group is shrunk
#' towards 0 by the method of moments: `tau2` is the variance between the
#' group means beyond what the noise explains.
#'
#' @param y The outcome: finite numbers, at least two cases. A `clpre` list
#'   from [clprep()] is taken too (the `taph-*` presets): its elements fill
#'   the arguments by name.
#' @param x The predictors: a numeric matrix or data frame of finite
#'   numbers, one row per case (code a factor as indicator columns first).
#' @param grp `NULL`, or the group of each case: labels, at least two
#'   groups.
#' @param pm,psd The prior means and standard deviations of the intercept
#'   and the coefficients: one number each or one per term (`psd`
#'   positive).
#' @inheritParams cntreg
#' @return A list: `coef`, a data frame with one row per term (`term`,
#'   `mean`, `sd`, the interval `lo`, `hi`, `ppos` the posterior
#'   probability that the coefficient is above 0, `prior` the prior's
#'   standard deviation); `sigma`, the noise standard deviation; `tau2`;
#'   `groups`, `NULL` or a data frame with one row per group (`group`, `n`,
#'   `raw` the mean residual, `shrink` the factor, `pooled` the shrunk
#'   intercept); `fitted`; `level`.
#' @references
#'   Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A.
#'   and Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd ed. CRC Press.
#' @examples
#' set.seed(1)
#' d <- data.frame(lime = rbinom(30, 1, 0.5), temp = rnorm(30, 12, 3))
#' y <- 0.3 + 0.2 * d$lime - 0.01 * d$temp + rnorm(30, sd = 0.05)
#' tpbhm(y, d, grp = rep(c("a", "b", "c"), 10), pm = c(0, 0.2, 0))$coef
#' @export
tpbhm <- function(y, x, grp = NULL, pm = 0, psd = 10, level = 0.95) {
  if (inherits(y, "clpre")) {
    return(.prep(tpbhm, y, match.call(), parent.frame()))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_tpbhm, y, x, grp, pm, psd, level)
}
