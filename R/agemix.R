#' The age-crime curve of a mixture of offender types
#'
#' The aggregate age-crime curve is the share-weighted mean of the curves
#' of latent types (adolescence-limited, life-course persistent, and so
#' on), and one aggregate curve fits many mixtures: the aggregate alone
#' cannot tell one type from several (Moffitt 1993; Gottfredson and
#' Hirschi 1990). What does tell them apart is the covariance of an
#' individual's offending at two ages: with one type it is 0 (given the
#' type, ages are independent), with several it is the covariance across
#' types of their curves at the two ages, which `cov` gives for every pair
#' of ages. With a panel of individuals observed at two ages, the observed
#' covariance and its permutation test say whether one type is enough.
#' See `proof("agemix")`.
#'
#' @param curves `NULL`, or the curves of the types: a numeric matrix or a
#'   data frame, one row per age and one column per type (at most 30).
#' @param shares The share of each type: numbers from 0 that sum to 1.
#' @param ages `NULL` (1, 2, ...) or the age of each row of `curves`.
#' @param panel `NULL`, or the outcomes of individuals at two ages: a
#'   numeric matrix or data frame of two columns, at least three rows.
#' @param nperm The number of permutations of the second age's outcomes for
#'   the panel test: a whole number from 0.
#' @param seed `NULL` (the draws follow [set.seed()]) or a whole number
#'   from 0.
#' @return An object of class `agemix`, a list: `curve` (`age`, `agg` the
#'   aggregate, then the types' curves; with `curves`), `cov` (the
#'   covariance across types between ages, a matrix named by age; with
#'   `curves`) and `panel` (`cov`, `p` the permutation p-value, `NA`
#'   without permutations, `n`, `nperm`; with `panel`); the parts not asked
#'   for are `NULL`.
#' @references
#'   Moffitt, T. E. (1993). Adolescence-limited and life-course-persistent
#'   antisocial behavior: a developmental taxonomy. *Psychological Review*,
#'   100(4), 674-701. \doi{10.1037/0033-295X.100.4.674}
#' @examples
#' ages <- 12:30
#' al <- dnorm(ages, 17, 2)
#' lcp <- rep(0.08, length(ages))
#' r <- agemix(cbind(al = al, lcp = lcp), c(0.8, 0.2), ages)
#' r$cov["15", "25"]
#' @export
agemix <- function(
  curves = NULL,
  shares = NULL,
  ages = NULL,
  panel = NULL,
  nperm = 2000,
  seed = NULL
) {
  r <- .Call(C_clab_r_agemix, curves, shares, ages, panel, nperm, seed)
  structure(r, class = "agemix")
}

#' Cumulative risk of a first event by age
#'
#' The share of a birth cohort that has had a first arrest, conviction or
#' prison admission by each age, from first events by age and the
#' population: the hazard at each age is the first events over those
#' still without one, and the cumulative risk is one minus the product of
#' the complements, with Greenwood's interval on the log scale of the
#' survival (Bonczar 2003; Pettit and Western 2004). With death rates the
#' first events count only among the cohort still alive, and no interval
#' is given.
#'
#' @param age The ages, increasing.
#' @param ev The first events at each age: numbers from 0.
#' @param pop The population at each age: positive numbers.
#' @param atrisk `NULL` (the population less those already with an event)
#'   or the number at risk at each age.
#' @param dth `NULL`, or the death rate at each age, from 0 to below 1.
#' @param level The level of the interval: one number between 0 and 1.
#' @return A data frame of class `cumrsk`: `age`, `hazard`, `atrisk`,
#'   `cumrisk`, `lo`, `hi` (`NA` with death rates).
#' @references Pettit, B. and Western, B. (2004). Mass imprisonment and the
#'   life course: race and class inequality in U.S. incarceration.
#'   *American Sociological Review*, 69(2), 151-169.
#'   \doi{10.1177/000312240406900201}
#' @examples
#' cumrsk(18:22, c(120, 100, 80, 60, 40), rep(10000, 5))
#' @export
cumrsk <- function(age, ev, pop, atrisk = NULL, dth = NULL, level = 0.95) {
  r <- .Call(C_clab_r_cumrsk, age, ev, pop, atrisk, dth, level)
  structure(r, class = c("cumrsk", "data.frame"))
}

#' Rates between stages of the justice pipeline
#'
#' Cases moving from one stage to the next (charges to convictions,
#' admissions to releases) take time, so a rate of this period's outflow
#' over this period's inflow mixes cohorts. `jpipe()` divides the
#' numerator summed over a window by the denominator summed over the same
#' window moved back by the lag between the stages.
#'
#' @param num The counts of the later stage, in time order: numbers from
#'   0.
#' @param den The counts of the earlier stage: numbers from 0, one per
#'   `num`.
#' @param win The window, in periods: a whole number from 1.
#' @param lag The lag between the stages, in periods: a whole number from
#'   0.
#' @param grp `NULL`, or the group of each period (each group in time
#'   order): labels.
#' @return A data frame of class `jpipe` with one row per period, the
#'   groups in label order: `group`, `index` (the period within its group),
#'   `num` and `den` (the window sums, `NA` until a full window), `rate`.
#' @examples
#' jpipe(num = c(5, 7, 6, 9, 8, 10), den = c(20, 22, 19, 25, 24, 26),
#'       win = 3, lag = 1)
#' @export
jpipe <- function(num, den, win = 12, lag = 0, grp = NULL) {
  r <- .Call(C_clab_r_jpipe, num, den, win, lag, grp)
  structure(r, class = c("jpipe", "data.frame"))
}

#' @export
print.agemix <- function(x, ...) {
  cat("agemix: age-crime curves of a mixture of types\n")
  if (!is.null(x[["curve"]])) {
    print(x[["curve"]], row.names = FALSE)
  }
  if (!is.null(x[["panel"]])) {
    print(x[["panel"]])
  }
  invisible(x)
}

#' @export
print.cumrsk <- function(x, ...) {
  cat("cumrsk: cumulative risk of a first event by age\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

#' @export
print.jpipe <- function(x, ...) {
  cat("jpipe: lagged rates between stages\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}
