#' Fairness metrics of a binary risk tool by group
#'
#' The confusion matrix of a risk score cut at `cut` within each group, and
#' the usual rates computed from it: accuracy, sensitivity (true positive
#' rate), specificity, false positive and false negative rates, positive and
#' negative predictive values, the share and number predicted positive and
#' the Matthews correlation, each with its ratio to the base group and, for
#' the first eight, a Wilson score interval (Wilson 1927). A case is
#' predicted positive when its score exceeds `cut`. No single metric
#' defines fairness, and when base rates differ between groups the error
#' rates and the predictive values cannot all be equal at once
#' (Chouldechova 2017); the table shows which differ, not which should not.
#'
#' @param score Risk scores: finite numbers, one per case.
#' @param y Outcomes: 0 or 1 (or `FALSE`/`TRUE`), one per case.
#' @param grp Group labels: numbers, text or a factor, one per case, no
#'   missing values. Groups are ordered as the factor's levels, or sorted.
#' @param cut The cut: one number.
#' @param base The group the ratios are taken against: its label, or its
#'   position from 1.
#' @param conf The confidence level of the intervals, strictly between 0
#'   and 1.
#' @return A list of data frames, one row per group with its label in
#'   `group`: `cnt` (`tp`, `fp`, `tn`, `fn`, `n`), `met` (`acc`, `tpr`,
#'   `tnr`, `fpr`, `fnr`, `ppv`, `npv`, `ppos`, `npos`, `mcc`; `NaN` where a
#'   rate has an empty denominator), `rat` (each metric over the base
#'   group's), `lo` and `hi` (the interval ends of the first eight metrics;
#'   the last two columns `NaN`).
#' @references Wilson, E. B. (1927). Probable inference, the law of
#'   succession, and statistical inference. *Journal of the American
#'   Statistical Association*, 22(158), 209-212.
#'   \doi{10.1080/01621459.1927.10502953}
#'
#'   Chouldechova, A. (2017). Fair prediction with disparate impact: a study
#'   of bias in recidivism prediction instruments. *Big Data*, 5(2),
#'   153-163. \doi{10.1089/big.2016.0047}
#' @examples
#' set.seed(1)
#' g <- rep(c("a", "b"), each = 100)
#' y <- rbinom(200, 1, ifelse(g == "a", 0.3, 0.5))
#' s <- plogis(rnorm(200, ifelse(y == 1, 0.8, -0.8)))
#' fair(s, y, g)$met
#' @export
fair <- function(score, y, grp, cut = 0.5, base = 1, conf = 0.95) {
  .Call(C_clab_r_fair, score, y, grp, cut, base, conf)
}

#' Two-level model for multivariate trace evidence
#'
#' `tlvc()` estimates the background of the two-level normal model of
#' Aitken and Lucy (2004) from a reference collection, `x_ij = theta_i +
#' e_ij`, items `theta_i ~ N(mu, C)` and measurements within an item
#' `e_ij ~ N(0, U)`: the overall mean, the pooled within-item covariance
#' and the ANOVA between-item covariance (unbiased for any pattern of
#' replication). `tlvlr()` gives the likelihood ratio of same source against
#' different sources for two items (a control and a recovered fragment, say)
#' from their mean vectors, computed in logarithms from Cholesky factors so
#' that atypical or distant items give a finite value.
#'
#' @param x The reference measurements: a numeric matrix, one row per
#'   observation and one column per variable.
#' @param item The item of each row: labels, every item at least once, some
#'   item more than once.
#' @param xa,xb The mean vectors of the two items compared: one number per
#'   variable.
#' @param na,nb Their numbers of replicates: positive.
#' @param bg The background, as `tlvc()` returns it.
#' @return `tlvc()`: a list of `mu` (the mean), `u` and `c` (the within- and
#'   between-item covariance matrices; `c` need not be positive definite,
#'   which `tlvlr()` then refuses). `tlvlr()`: a named vector, `ln` and
#'   `log10`, the natural and base-10 logarithms of the ratio.
#' @references Aitken, C. G. G. and Lucy, D. (2004). Evaluation of trace
#'   evidence in the form of multivariate data. *Journal of the Royal
#'   Statistical Society, Series C*, 53(1), 109-122.
#'   \doi{10.1046/j.0035-9254.2003.05271.x}
#' @examples
#' set.seed(1)
#' th <- matrix(rnorm(20, 10, 2), 10, 2)
#' x <- th[rep(1:10, each = 3), ] + matrix(rnorm(60, 0, 0.2), 30, 2)
#' bg <- tlvc(x, rep(1:10, each = 3))
#' tlvlr(th[1, ], th[1, ] + 0.05, 3, 3, bg)
#' tlvlr(th[1, ], th[2, ], 3, 3, bg)
#' @export
tlvc <- function(x, item) {
  .Call(C_clab_r_tlvc, x, item)
}

#' @rdname tlvc
#' @export
tlvlr <- function(xa, xb, na, nb, bg) {
  .Call(C_clab_r_tlvlr, xa, xb, na, nb, bg)
}

#' Aoristic weights
#'
#' Spreads each event known only to have happened in an interval `[from,
#' to)` with uniform density over that interval and counts each instant in
#' the bin of its local clock time (Ratcliffe 2002): with the defaults, the
#' 168 hours of the week from Sunday 00:00. A change of clock (daylight
#' saving) moves the remaining time to the right hours: skipped hours get
#' nothing and a repeated hour gets both. Missing values follow stated
#' rules: an event with no `from` is left out; with no `to`, a reversed
#' interval or one shorter than `lo` seconds, its whole weight goes to the
#' start's bin; one of at least `hi` seconds is spread evenly over all bins.
#'
#' @param from,to The start and end of each event's interval: date-times
#'   (POSIXct) or seconds since 1970-01-01 UTC; `NA` allowed as above, not
#'   infinite values.
#' @param bw The bin width in seconds, positive.
#' @param nb The number of bins in a period: one whole number from 1.
#' @param origin The local clock time at which bin 1 starts, in the same
#'   units (the default is Sunday 1970-01-04 00:00).
#' @param off The offset of the local clock from UTC in seconds, before
#'   `tt[1]` (for a fixed offset, the only one).
#' @param tt,toff `NULL`, or the instants at which the local clock's
#'   offset changes (strictly increasing) and the new offsets.
#' @param lo,hi The lengths in seconds below which an interval counts at its
#'   start and from which it counts everywhere, `0 <= lo < hi` (`lo = 0`,
#'   `hi = Inf` spread every interval exactly).
#' @param each `TRUE` to return each event's weights.
#' @return A list: `tot` (the summed weights by bin), `w` (`NULL`, or a
#'   matrix with one row per event and one column per bin, rows summing to
#'   1, 0 for events left out) and `cnt` (events left out, without end,
#'   reversed, short and long).
#' @references Ratcliffe, J. H. (2002). Aoristic signatures and the
#'   spatio-temporal analysis of high volume crime patterns. *Journal of
#'   Quantitative Criminology*, 18(1), 23-43.
#'   \doi{10.1023/A:1013240828824}
#' @examples
#' f <- as.POSIXct(c("2024-03-04 22:00", "2024-03-05 09:30"), tz = "UTC")
#' t <- as.POSIXct(c("2024-03-05 02:00", "2024-03-05 10:30"), tz = "UTC")
#' a <- aorwt(f, t, bw = 3600, nb = 24, origin = 0)
#' a$tot
#' @export
aorwt <- function(from, to, bw = 3600, nb = 168, origin = 259200, off = 0,
                  tt = NULL, toff = NULL, lo = 120, hi = 604800,
                  each = FALSE) {
  .Call(C_clab_r_aorwt, from, to, bw, nb, origin, off, tt, toff, lo, hi,
        each)
}

#' Benford digit screens
#'
#' `bdig()` takes the first `nd` significant digits of each value. With
#' `dp = -1` they are the digits of the shortest decimal that gives back the
#' double, decided exactly, so 0.29 gives 29 where truncating `x` times a
#' power of ten gives 28; with `dp` from 0 to 22 the data are taken as
#' recorded to `dp` decimals. `bnft()` tests the digits against Benford's
#' law, `P(d) = log10(1 + 1 / d)`, by the chi-square statistic with its
#' exact upper tail, a Monte Carlo p-value from samples of Benford digits of
#' the same size, and Nigrini's mean absolute deviation. A departure from
#' Benford's law is a reason to look, not a finding of fabrication: many
#' honest data sets (bounded amounts, assigned numbers) do not follow it.
#'
#' @param x Numbers. Zeros, non-finite values (including `NA`), values of
#'   the excluded sign and values outside the exact range are set aside and
#'   counted.
#' @param nd The number of leading digits: from 1 to 9 for `bdig()`, 1 to
#'   6 for `bnft()`.
#' @param dp `-1`, or the number of decimals the data were recorded to.
#' @param sgn Which values to keep: `"pos"`, `"neg"` or `"both"` (by
#'   absolute value).
#' @param dig Digits as `bdig()` returns them (0 is skipped).
#' @param nsim The number of Monte Carlo samples: a whole number from 0.
#' @inheritParams nnidx
#' @return `bdig()`: a list of `dig` (an integer per value, 0 when set
#'   aside) and `cnt` (zeros, non-finite, excluded sign, out of range).
#'   `bnft()`: a list of `obs` (a data frame of `digit`, observed `count` and
#'   Benford `exp`ected count) and `stat` (`m` digits used, `chisq`, `df`,
#'   `p`, `pmc` (Monte Carlo, `NaN` without samples) and `mad`).
#' @references Benford, F. (1938). The law of anomalous numbers.
#'   *Proceedings of the American Philosophical Society*, 78(4), 551-572.
#'
#'   Nigrini, M. J. (2012). *Benford's Law: Applications for Forensic
#'   Accounting, Auditing, and Fraud Detection*. Wiley.
#'   \doi{10.1002/9781119203094}
#' @examples
#' set.seed(1)
#' x <- exp(runif(500, 0, 10))
#' d <- bdig(x)$dig
#' bnft(d, nsim = 199, seed = 1)$stat
#' @export
bdig <- function(x, nd = 1, dp = -1, sgn = "pos") {
  .Call(C_clab_r_bdig, x, nd, dp, sgn)
}

#' @rdname bdig
#' @export
bnft <- function(dig, nd = 1, nsim = 0, seed = NULL) {
  .Call(C_clab_r_bnft, dig, nd, nsim, seed)
}

#' Scott-Knott grouping of measurements
#'
#' Splits sorted measurements into groups that differ (Scott and Knott
#' 1974) when the variance of one measurement is known, as in the forensic
#' glass version of Triggs, Curran, Buckleton and Walsh (1997): a segment is
#' split at the cut with the largest between-groups sum of squares when that
#' is significant at `alpha`, judged against `nsim` simulated samples of
#' independent normal values with the known variance (a test of exact size
#' under that null).
#'
#' @param x Measurements: finite numbers (refractive indices, say).
#' @param s2 The variance of one measurement, on the scale of `x`: positive.
#' @param alpha The level of each split's test, strictly between 0 and 1.
#' @param nsim The number of simulated samples per test: a whole number from
#'   1.
#' @inheritParams nnidx
#' @return An integer vector of group labels, one per measurement in input
#'   order, numbered from 1 in increasing order of the values.
#' @references Scott, A. J. and Knott, M. (1974). A cluster analysis method
#'   for grouping means in the analysis of variance. *Biometrics*, 30(3),
#'   507-512. \doi{10.2307/2529204}
#'
#'   Triggs, C. M., Curran, J. M., Buckleton, J. S. and Walsh, K. A. J.
#'   (1997). The grouping problem in forensic glass analysis: a divisive
#'   approach. *Forensic Science International*, 85(1), 1-14.
#'   \doi{10.1016/S0379-0738(96)02037-3}
#' @examples
#' set.seed(1)
#' x <- c(rnorm(5, 1.5180, 4e-5), rnorm(5, 1.5195, 4e-5))
#' skp(x, s2 = 4e-5^2, nsim = 199, seed = 1)
#' @export
skp <- function(x, s2, alpha = 0.05, nsim = 999, seed = NULL) {
  .Call(C_clab_r_skp, x, s2, alpha, nsim, seed)
}

#' Constant reporters of a reporting panel
#'
#' The units (agencies reporting to a national crime statistics programme,
#' say) that report at least `minrep` months in every year from `y0` to
#' `y1`. A year without a record of the unit, or with too few months,
#' disqualifies it, whatever other units do that year. Counts over time are
#' comparable only on such a constant panel; a trend from all reporters
#' mixes change in crime with change in who reports.
#'
#' @param unit The unit of each record: labels, no missing values.
#' @param year The year of each record: whole numbers.
#' @param rep The months reported in each record: whole numbers.
#' @param y0,y1 The first and last years, `y0 <= y1`, from 0.
#' @param minrep The fewest months that count as reporting: a whole number
#'   from 0.
#' @return A data frame with one row per unit (sorted): `unit` and `const`
#'   (`TRUE` for a constant reporter).
#' @references Maltz, M. D. (1999). *Bridging Gaps in Police Crime Data*.
#'   Bureau of Justice Statistics, NCJ 176365.
#'   <https://bjs.ojp.gov/content/pub/pdf/bgpcd.pdf>
#' @examples
#' rptc(unit = c("a", "a", "b", "b", "b"), year = c(2020, 2021, 2020, 2021,
#'      2021), rep = c(12, 12, 12, 6, 12), y0 = 2020, y1 = 2021)
#' @export
rptc <- function(unit, year, rep, y0, y1, minrep = 12) {
  .Call(C_clab_r_rptc, unit, year, rep, y0, y1, minrep)
}

#' Rosenbaum sensitivity bounds for matched pairs
#'
#' How strong would hidden bias have to be to explain away a matched-pairs
#' result? `gam` is Gamma, the most the odds of treatment can differ within
#' a pair because of an unobserved covariate (1: no hidden bias).
#' `rswp()` gives the bounds on the one-sided p-value of Wilcoxon's
#' signed-rank test at each Gamma (Rosenbaum 1987), from the normal
#' approximation, and the Gamma at which the upper bound reaches `alpha`.
#' `rshl()` gives the bounds on the Hodges-Lehmann estimate of the effect
#' (Rosenbaum 1993), computed exactly from the Walsh averages, with no
#' grid.
#'
#' @param d The treated-minus-control differences: finite numbers, one per
#'   pair.
#' @param gam Values of Gamma: numbers from 1.
#' @param alpha The level for the critical Gamma, strictly between 0 and
#'   0.5.
#' @return `rswp()`: a list of `p` (a data frame of `gam`, `plo`, `pup`),
#'   `gcrit` (the critical Gamma, below 1 when the test is not significant
#'   without bias) and `t` (the signed-rank statistic). `rshl()`: a data
#'   frame of `gam`, `lo`, `hi`. `NaN` when every difference is 0.
#' @references Rosenbaum, P. R. (1987). Sensitivity analysis for certain
#'   permutation inferences in matched observational studies. *Biometrika*,
#'   74(1), 13-26. \doi{10.1093/biomet/74.1.13}
#'
#'   Rosenbaum, P. R. (1993). Hodges-Lehmann point estimates of treatment
#'   effect in observational studies. *Journal of the American Statistical
#'   Association*, 88(424), 1250-1253. \doi{10.1080/01621459.1993.10476405}
#' @examples
#' set.seed(1)
#' d <- rnorm(50, 0.5)
#' rswp(d, gam = c(1, 1.5, 2))
#' rshl(d, gam = c(1, 1.5, 2))
#' @export
rswp <- function(d, gam = 1, alpha = 0.05) {
  .Call(C_clab_r_rswp, d, gam, alpha)
}

#' @rdname rswp
#' @export
rshl <- function(d, gam = 1) {
  .Call(C_clab_r_rshl, d, gam)
}

#' Group-based trajectory model for counts
#'
#' The finite mixture of Nagin (1999) for counts over time, zero-inflated
#' Poisson within groups (Lambert 1992; Nielsen et al. 2014): unit `i`
#' belongs to group `k` with probability `prob[k]`; given `k`, its counts
#' are independent Poisson with log mean `off + x_t' beta_k` and, when `z` is
#' given, zero with extra probability `plogis(z_t' gamma_k)`. `zgll()`
#' evaluates the log-likelihood and the posterior memberships; `zgem()`
#' fits by a generalised EM from the values given (the log-likelihood never
#' decreases). One start per call: the model has local maxima, so several
#' starts are the caller's loop, and the number of groups is the caller's
#' choice (compare the log-likelihoods with a penalty such as BIC).
#'
#' @param y Counts: a matrix with one row per unit and one column per
#'   period; `NA` or a negative value is a missing period.
#' @param x The design of the Poisson part: a matrix with one row per
#'   period (columns 1, t, t^2, say).
#' @param z `NULL` for a Poisson model, or the design of the zero-inflation
#'   part, one row per period.
#' @param par The coefficients: a matrix with one row per column of `x`
#'   (then of `z`) and one column per group.
#' @param prob The group shares: positive numbers summing to 1.
#' @param off `NULL`, or log exposures with the dimensions of `y`.
#' @param maxit The most iterations: a whole number from 1.
#' @param tol Stop when an iteration gains at most `tol (1 + |ll|)`.
#' @return `zgll()`: a list of `ll` and `post` (a matrix of posterior
#'   memberships, one row per unit). `zgem()`: a list of the fitted `par`,
#'   `prob`, `post`, `ll` and `iter` (iterations made).
#' @references Nagin, D. S. (1999). Analyzing developmental trajectories: a
#'   semiparametric, group-based approach. *Psychological Methods*, 4(2),
#'   139-157. \doi{10.1037/1082-989X.4.2.139}
#'
#'   Lambert, D. (1992). Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
#'   \doi{10.2307/1269547}
#' @examples
#' set.seed(1)
#' x <- cbind(1, (1:6 - 3.5) / 2.5)
#' b <- cbind(c(0.5, -1), c(1.5, 0.5))
#' y <- t(sapply(rep(1:2, each = 50),
#'               function(k) rpois(6, exp(x %*% b[, k]))))
#' f <- zgem(y, x, par = b + 0.2, prob = c(0.5, 0.5))
#' f$par
#' @export
zgll <- function(y, x, z = NULL, par, prob, off = NULL) {
  .Call(C_clab_r_zgll, y, x, z, par, prob, off)
}

#' @rdname zgll
#' @export
zgem <- function(y, x, z = NULL, par, prob, off = NULL, maxit = 500,
                 tol = 1e-8) {
  .Call(C_clab_r_zgem, y, x, z, par, prob, off, maxit, tol)
}

#' Group-based multivariate trajectories
#'
#' The group-based trajectory model with several normal indicators per unit
#' (Nagin 1999; Magrini 2022): given group `k`, the `r` indicators of a unit
#' at period `t` are multivariate normal with mean `B_k' x_t` and covariance
#' `S_k`, independently over periods. `mvtll()` evaluates the
#' log-likelihood and posterior memberships; `mvtem()` fits by EM with
#' closed-form steps from the values given. A period with any missing
#' indicator is left out for that unit.
#'
#' @param y The indicators: an array of units by periods by indicators;
#'   `NA` allowed, not infinite values.
#' @param x The design: a matrix with one row per period.
#' @param beta The coefficients: an array of design columns by indicators by
#'   groups.
#' @param sig The covariances: an array of indicators by indicators by
#'   groups, each positive definite.
#' @inheritParams zgll
#' @return `mvtll()`: a list of `ll` and `post`. `mvtem()`: a list of the
#'   fitted `beta`, `sig`, `prob`, `post`, `ll` and `iter`.
#' @references Magrini, A. (2022). Assessment of agricultural sustainability
#'   in European Union countries: a group-based multivariate trajectory
#'   approach. *AStA Advances in Statistical Analysis*, 106.
#'   \doi{10.1007/s10182-022-00437-9}
#' @examples
#' set.seed(1)
#' x <- cbind(1, 1:5)
#' y <- array(rnorm(40 * 5 * 2), c(40, 5, 2))
#' y[1:20, , 1] <- y[1:20, , 1] + 3
#' b <- array(c(0, 0, 0, 0, 3, 0, 0, 0), c(2, 2, 2))
#' s <- array(diag(2), c(2, 2, 2))
#' mvtem(y, x, b, s, prob = c(0.5, 0.5))$prob
#' @export
mvtll <- function(y, x, beta, sig, prob) {
  .Call(C_clab_r_mvtll, y, x, beta, sig, prob)
}

#' @rdname mvtll
#' @export
mvtem <- function(y, x, beta, sig, prob, maxit = 500, tol = 1e-8) {
  .Call(C_clab_r_mvtem, y, x, beta, sig, prob, maxit, tol)
}

#' Glass fragments recovered under transfer and persistence
#'
#' The probability of recovering 0, 1, ... fragments of glass from a
#' suspect under the transfer, persistence and recovery model of Curran,
#' Triggs, Buckleton, Walsh and Hicks (1998), with uncertain parameters
#' averaged over: the transfer amount and the distance deterministically
#' (by quadrature), the time since the offence exactly, and the loss and
#' recovery rates by `nsim` Monte Carlo draws (exactly when every rate is
#' fixed).
#'
#' @param par The 17 parameters, in order: `lam`, `lsd` (transfer amount
#'   N(lam, lsd^2), used as its absolute value), `d` (distance scale), `Q`
#'   (share of high-persistence fragments), `l0`, `u0` (first-hour loss of
#'   the others, uniform), `ls0`, `us0` (of the high-persistence ones),
#'   `lj`, `uj` and `lsj`, `usj` (each later hour), `lR`, `uR` (share
#'   recovered), `lt`, `ut` (hours since the offence: fixed when equal) and
#'   `r` (the negative binomial probability of the time).
#' @param deffect `TRUE` for the distance effect, `FALSE` without it.
#' @param tdist The law of the time when `lt < ut`: `"nb"` (negative
#'   binomial with size `(lt + ut) / 2` and probability `r`), `"nbtrunc"`
#'   (that law on `floor(lt)` to `ceiling(ut)`) or `"unif"` (uniform on those
#'   hours).
#' @param ymax The largest count reported: a whole number from 0.
#' @param nsim The Monte Carlo draws of the rates: a whole number from 1.
#' @inheritParams nnidx
#' @param hl,nq The quadrature's grid step (1e-4 to 1) and panels (1 to
#'   100000); the defaults resolve the model's usual settings.
#' @return A list: `prob`, a data frame of `y` (0 to `ymax`), `p` and `se`
#'   (the Monte Carlo standard error, 0 when exact); `cut`, the time mass
#'   left out and the intensity mass on the grid.
#' @references Curran, J. M., Triggs, C. M., Buckleton, J. S., Walsh, K. A.
#'   J. and Hicks, T. (1998). Assessing transfer probabilities in a Bayesian
#'   interpretation of forensic glass evidence. *Science and Justice*,
#'   38(1), 15-21. \doi{10.1016/S1355-0306(98)72068-4}
#' @examples
#' par <- c(80, 50, 0.5, 0.05, 0.8, 0.9, 0.1, 0.2, 0.4, 0.6, 0.05, 0.1,
#'          0.5, 0.7, 1, 1, 0.5)
#' tfrp(par, ymax = 10, nsim = 100, seed = 1)$prob
#' @export
tfrp <- function(par, deffect = TRUE, tdist = "nb", ymax = 20, nsim = 1000,
                 seed = NULL, hl = 0.02, nq = 256) {
  .Call(C_clab_r_tfrp, par, deffect, tdist, ymax, nsim, seed, hl, nq)
}

#' Lowess smoothing
#'
#' Robust locally weighted regression (Cleveland 1979) with the conventions
#' of the classic implementation, so that the fit equals that of
#' [stats::lowess()]; the fitted values come back in the order of the
#' input, not sorted.
#'
#' @param x,y The points: finite numbers of equal length.
#' @param f The share of points in each neighbourhood: positive.
#' @param iter The number of robustness steps: a whole number from 0.
#' @param delta Points within `delta` of the last fitted one are
#'   interpolated: `NULL` for 1 per cent of the range of `x`, or a number
#'   from 0.
#' @return A numeric vector of fitted values, one per point in input order.
#' @references Cleveland, W. S. (1979). Robust locally weighted regression
#'   and smoothing scatterplots. *Journal of the American Statistical
#'   Association*, 74(368), 829-836. \doi{10.1080/01621459.1979.10481038}
#' @examples
#' lwss(c(3, 1, 2, 5, 4), c(2, 1, 3, 6, 4))
#' @export
lwss <- function(x, y, f = 2 / 3, iter = 3, delta = NULL) {
  .Call(C_clab_r_lwss, x, y, f, iter, delta)
}

#' The Chumbley score for two tool-mark profiles
#'
#' Compares two striated tool-mark profiles by the deterministic Chumbley
#' score (Chumbley et al. 2010; Hadler and Morris 2018): each profile is
#' trimmed and detrended by its lowess smooth, the pair of windows of
#' `wopt` points with the largest correlation fixes the shift, and windows
#' of `wval` points at that shift are compared with windows at other
#' shifts by a standardised Wilcoxon rank-sum statistic.
#'
#' @param a,b The two profiles: finite numbers (heights along the mark).
#' @param wopt,wval The window lengths of the optimisation and validation
#'   steps: whole numbers from 2 with `wval <= wopt`.
#' @param coarse The lowess span of the detrending: positive.
#' @return A named numeric vector: `nsame`, `ndiff` (window pairs used),
#'   `u` (the statistic) and `p` (its upper normal tail; `NaN` when either
#'   set is empty).
#' @references Chumbley, L. S., Morris, M. D., Kreiser, M. J. et al.
#'   (2010). Validation of tool mark comparisons obtained using a quantitative,
#'   comparative, statistical algorithm. *Journal of Forensic Sciences*,
#'   55(4), 953-961. \doi{10.1111/j.1556-4029.2010.01424.x}
#'
#'   Hadler, J. R. and Morris, M. D. (2018). An improved version of a
#'   tool mark comparison algorithm. *Journal of Forensic Sciences*, 63(3),
#'   849-855. \doi{10.1111/1556-4029.13640}
#' @examples
#' set.seed(1)
#' m <- cumsum(rnorm(700))
#' a <- m[1:600] + rnorm(600, sd = 0.2)
#' b <- m[51:650] + rnorm(600, sd = 0.2)
#' chmb(a, b, wopt = 120, wval = 30)
#' @export
chmb <- function(a, b, wopt, wval, coarse = 0.25) {
  .Call(C_clab_r_chmb, a, b, wopt, wval, coarse)
}

#' Local Getis-Ord Gi* on a grid of counts
#'
#' The local statistic Gi* of Getis and Ord (1992; Ord and Getis 1995) for
#' each cell of a grid of counts, with binary weights on the cells within
#' `rad` cell widths (`rad = 1.5`: the 3 by 3 block, cell included), its
#' normal p-value and Benjamini-Hochberg q-value, a conditional permutation
#' p-value, and an exposure-aware version: given the total, the count of a
#' neighbourhood is binomial with the neighbourhood's share of the exposure
#' under complete spatial randomness, so cells overhanging the study area
#' neither inflate the mean nor get flagged. No threshold is applied.
#'
#' @param cnt Counts per cell: a matrix of whole numbers from 0.
#' @param ex `NULL` for equal exposure, or a matrix of the same dimensions:
#'   each cell's exposure (its area inside the study window, say), 0 for a
#'   cell outside (its results are `NaN`).
#' @param rad The neighbourhood radius in cell widths, from 0 to 1000 (1:
#'   the cell and its four rook neighbours; 1.5: the queen neighbourhood).
#' @param nsim The number of permutations: a whole number from 0.
#' @inheritParams nnidx
#' @return A data frame with one row per cell (down the columns of `cnt`):
#'   `row`, `col`, `gz` (the classic Gi* z), `pz`, `qz`, `pp` (the
#'   permutation p-value, `NaN` without permutations), `ze` (the
#'   exposure-aware z), `pe` (its exact two-sided binomial p-value) and
#'   `qe`.
#' @references Getis, A. and Ord, J. K. (1992). The analysis of spatial
#'   association by use of distance statistics. *Geographical Analysis*,
#'   24(3), 189-206. \doi{10.1111/j.1538-4632.1992.tb00261.x}
#'
#'   Ord, J. K. and Getis, A. (1995). Local spatial autocorrelation
#'   statistics: distributional issues and an application. *Geographical
#'   Analysis*, 27(4), 286-306. \doi{10.1111/j.1538-4632.1995.tb00912.x}
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(100, 2), 10)
#' m[4:5, 4:5] <- 9
#' g <- gstr(m)
#' g[g$qe < 0.05, ]
#' @export
gstr <- function(cnt, ex = NULL, rad = 1.5, nsim = 0, seed = NULL) {
  .Call(C_clab_r_gstr, cnt, ex, rad, nsim, seed)
}

#' Andresen's spatial point pattern test
#'
#' Compares where two point patterns fall over the same areal units
#' (Andresen 2009): for each unit, the difference between its share of the
#' test pattern `b` and its share of the base pattern `a`, with a bootstrap
#' p-value from multinomial resamples of each pattern; the S index is the
#' share of units whose difference is not significant (1: the patterns are
#' alike everywhere).
#'
#' @param a,b The counts of the base and the test pattern by unit: whole
#'   numbers from 0; `NA` in either leaves the unit out (counted, never
#'   taken as 0).
#' @param nboot The number of bootstrap resamples: a whole number from 1.
#' @inheritParams nnidx
#' @param alpha The level of each unit's test, strictly between 0 and 1.
#' @param scale The multiplier of the differences (100: percentage points).
#' @return A list: `units`, a data frame with one row per unit of `diff`, `p`
#'   and `flag` (the sign of a significant difference, 0 otherwise; `NaN`
#'   without a p-value); `glob`, a named vector: units `left` out, `empty` in
#'   both, `elig`ible and `flagged`, the S index `s`, its `z` against
#'   `1 - alpha`, the exact `p` of so many flags and `nboot`.
#' @references Andresen, M. A. (2009). Testing for similarity in area-based
#'   spatial patterns: a nonparametric Monte Carlo approach. *Applied
#'   Geography*, 29(3), 333-345. \doi{10.1016/j.apgeog.2008.12.004}
#' @examples
#' set.seed(1)
#' a <- rpois(30, 10)
#' b <- rpois(30, 10)
#' sptst(a, b, nboot = 199, seed = 1)$glob
#' @export
sptst <- function(a, b, nboot = 200, seed = NULL, alpha = 0.05,
                  scale = 100) {
  .Call(C_clab_r_sptst, a, b, nboot, seed, alpha, scale)
}
