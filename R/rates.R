#' Rates and proportions by group: ratios and differences
#'
#' The rate (events over exposure) or the proportion (events over trials)
#' of each group, and its ratio to and difference from a reference group,
#' with approximate and exact intervals. Rows that share a group are
#' pooled: give one row per group, or one per case (0/1 outcomes with
#' `n = 1`, the default for proportions), or anything between.
#'
#' For a rate `x / n` the approximate interval is `rate * exp(-+z /
#' sqrt(x))` and the exact one the Poisson interval through the gamma
#' quantiles, as [stats::poisson.test()]; the ratio's approximate interval
#' is log-normal with variance `1 / x + 1 / x0`, its exact interval the
#' one conditional on the total count (the Clopper-Pearson interval of `x
#' / (x + x0)`, as the two-sample [stats::poisson.test()]), and the
#' difference's is Wald. For a proportion the intervals are Wilson's
#' score interval and the exact Clopper-Pearson one ([stats::binom.test()]);
#' the ratio's is log-normal (Katz et al.), with no exact interval (`NA`),
#' and the difference's is Newcombe's hybrid score interval. With a zero
#' count the log-normal intervals are `NA`: no continuity correction is
#' made.
#'
#' `rep` is a group's share of the events over its share of the exposure
#' (or the trials): over 1 when the group is over-represented among the
#' events, as in a comparison of who is placed in segregation with who is
#' in custody.
#'
#' @param x The events: numbers from 0; or a `clpre` list from [clprep()].
#' @param n The exposure (person-time, population) of each row for rates,
#'   or the trials for proportions: positive numbers; `NULL` for one trial
#'   per row (proportions only).
#' @param grp `NULL` (each row is a group, named by the names of `x`) or
#'   the group of each row: labels.
#' @param type `"rate"` (events over exposure) or `"prop"` (a proportion).
#' @param per The rates and differences are given per `per` units of `n`
#'   (1e5 for rates per 100,000, 100 for percentages): one positive number.
#' @inheritParams cntreg
#' @return An object of class `grprr`, a list: `rates`, a data frame with
#'   one row per group (`group`, `x`, `n`, `rate`, `lo`, `hi` the
#'   approximate interval, `elo`, `ehi` the exact one, `rep`); `cmp`, one
#'   row per group but the reference (`group`, `ratio`, `lo`, `hi`, `elo`,
#'   `ehi`, `diff`, `dlo`, `dhi`); `ref`, `type`, `per` and `level`.
#' @references
#'   Katz, D., Baptista, J., Azen, S. P. and Pike, M. C. (1978). Obtaining
#'   confidence intervals for the risk ratio in cohort studies.
#'   *Biometrics*, 34(3), 469-474. \doi{10.2307/2530610}
#'
#'   Newcombe, R. G. (1998). Interval estimation for the difference between
#'   independent proportions: comparison of eleven methods. *Statistics in
#'   Medicine*, 17(8), 873-890.
#'   \doi{10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I}
#' @examples
#' # use-of-force incidents per 1,000 contacts by perceived race
#' grprr(c(white = 120, black = 95, other = 30),
#'       c(white = 60000, black = 21000, other = 15000), per = 1000)
#' # proportions from one row per case
#' grprr(c(1, 0, 0, 1, 1, 0, 1, 1), grp = c(1, 1, 1, 1, 2, 2, 2, 2),
#'       type = "prop")
#' @export
grprr <- function(
  x,
  n = NULL,
  grp = NULL,
  ref = NULL,
  type = "rate",
  per = 1,
  level = 0.95
) {
  if (inherits(x, "clpre")) {
    return(.prep(grprr, x, match.call(), parent.frame()))
  }
  .Call(C_clab_r_grprr, x, n, grp, ref, type, per, level)
}

#' Directly standardised rates
#'
#' The rate each group would have with the stratum make-up (age bands,
#' say) of a standard population: the sum over strata of the standard's
#' weights times the group's stratum rates. The interval is Fay and Feuer's
#' gamma interval; the ratio of two standardised rates has a log-normal
#' interval and their difference a normal one.
#'
#' @param strat The stratum of each row: labels.
#' @param std The standard population (or weights) of each stratum:
#'   numbers from 0 with a positive sum, named by stratum or in the order of
#'   the strata (a factor's levels, or the sorted labels).
#' @param grp `NULL` (one group) or the group of each row: labels. Each
#'   group has one row per stratum.
#' @param n The exposure of each row: positive numbers.
#' @inheritParams grprr
#' @return An object of class `dastd`, a list: `rates` (`group`, `x` the
#'   events, `rate`, `lo`, `hi`, `var` the variance of the rate), `cmp` (one
#'   row per group but the reference: `group`, `ratio`, `lo`, `hi`, `diff`,
#'   `dlo`, `dhi`), `ref`, `per` and `level`.
#' @references
#'   Fay, M. P. and Feuer, E. J. (1997). Confidence intervals for directly
#'   standardized rates: a method based on the gamma distribution.
#'   *Statistics in Medicine*, 16(7), 791-801.
#'   \doi{10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-#}
#' @examples
#' d <- data.frame(age = rep(c("18-29", "30-49", "50+"), 2),
#'                 grp = rep(c("A", "B"), each = 3),
#'                 x = c(30, 22, 5, 41, 30, 12),
#'                 n = c(5000, 9000, 6000, 3000, 8000, 7000))
#' dastd(d$x, d$n, d$age, d$grp, std = c(30, 45, 25), per = 1e4)
#' @export
dastd <- function(
  x,
  n,
  strat,
  grp = NULL,
  std,
  ref = NULL,
  per = 1,
  level = 0.95
) {
  if (inherits(x, "clpre")) {
    return(.prep(dastd, x, match.call(), parent.frame()))
  }
  .Call(C_clab_r_dastd, x, n, strat, grp, std, ref, per, level)
}

#' Indirect standardisation: standardised ratios with shrinkage
#'
#' The events observed in each group (an area, a police service) against
#' those expected from stratum rates, the standardised incidence ratio
#' (SIR) with its exact Poisson interval, funnel-plot limits for a ratio of
#' 1 at its expected count, and the ratio shrunk towards the overall one by
#' the empirical Bayes estimate of Marshall (a gamma prior fitted by
#' moments), so that groups with little exposure are not ranked on noise.
#'
#' @param rates `NULL` for the stratum rates pooled over all groups, or the
#'   reference rate of each stratum: numbers from 0, named by stratum or in
#'   the order of the strata.
#' @param grp The group of each row: labels.
#' @inheritParams dastd
#' @return An object of class `sirex`, a list: `sir`, a data frame (`group`,
#'   `obs`, `expd`, `sir`, `lo`, `hi`, `eb` the shrunk ratio, `shr` the
#'   share of the prior in it, `flo`, `fhi` the funnel limits, `out`
#'   `TRUE` outside them); `prior`, the prior's mean and variance (a
#'   variance of 0 when the ratios vary no more than chance, and then every
#'   `eb` is the mean); and `level`.
#' @references
#'   Marshall, R. J. (1991). Mapping disease and mortality rates using
#'   empirical Bayes estimators. *Journal of the Royal Statistical Society,
#'   Series C*, 40(2), 283-294. \doi{10.2307/2347593}
#' @examples
#' d <- data.frame(area = rep(c("N", "S", "E"), each = 2),
#'                 age = rep(c("young", "old"), 3),
#'                 x = c(12, 3, 30, 9, 4, 1),
#'                 n = c(4000, 2000, 8000, 3000, 1500, 900))
#' sirex(d$x, d$n, d$age, d$area)
#' @export
sirex <- function(x, n, strat, grp, rates = NULL, level = 0.95) {
  if (inherits(x, "clpre")) {
    return(.prep(sirex, x, match.call(), parent.frame()))
  }
  .Call(C_clab_r_sirex, x, n, strat, grp, rates, level)
}

#' @export
print.grprr <- function(x, ...) {
  cat(
    "grprr: ",
    if (x[["type"]] == "rate") "rates" else "proportions",
    if (x[["per"]] != 1) paste(" per", format(x[["per"]])),
    ", ",
    format(100 * x[["level"]]),
    "% intervals\n",
    sep = ""
  )
  print(x[["rates"]], row.names = FALSE)
  if (nrow(x[["cmp"]]) > 0L) {
    cat("against ", x[["ref"]], ":\n", sep = "")
    print(x[["cmp"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.dastd <- function(x, ...) {
  cat(
    "dastd: standardised rates",
    if (x[["per"]] != 1) paste(" per", format(x[["per"]])),
    "\n",
    sep = ""
  )
  print(x[["rates"]], row.names = FALSE)
  if (nrow(x[["cmp"]]) > 0L) {
    cat("against ", x[["ref"]], ":\n", sep = "")
    print(x[["cmp"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.sirex <- function(x, ...) {
  cat(
    "sirex: observed over expected, ",
    format(100 * x[["level"]]),
    "% intervals and funnel limits\n",
    sep = ""
  )
  print(x[["sir"]], row.names = FALSE)
  invisible(x)
}
