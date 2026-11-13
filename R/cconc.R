#' Concentration of crime across places, people or units
#'
#' How unequally a total is spread across units: the Gini coefficient and
#' the Lorenz curve, the share held by the top units, and what chance alone
#' would show. A few places holding most of the incidents, or a few people
#' most of the placements, is the law of crime concentration (Weisburd
#' 2015); with sparse counts much of the Gini comes from the units with no
#' events at all, and a Poisson process at the same mean already leaves
#' many units at zero (Bernasco and Steenbeek 2017).
#'
#' The Gini coefficient here is the population form, `1 - sum p (L[i-1] +
#' L[i])` over the units in ascending order (`p` a unit's share of the
#' units, `L` the cumulative share of the total), which for unit weights is
#' `2 sum i x[(i)] / (n T) - (n + 1) / n`. It splits exactly as `zero +
#' (1 - zero) gpos`, the share of units at zero plus the rest times the
#' Gini of the units with events; `gnull` is the same split at the share of
#' zeros a Poisson process at the observed mean leaves, `exp(-mean)`. The
#' variance over the mean (`disp`) is 1 for Poisson counts, and its excess
#' estimates the variance of the places' underlying rates (`ivar`). The
#' Gini of counts is at most `sqrt(vpop / 2) / mean` (`cvb`), so for
#' Poisson counts chance alone can produce a Gini up to `1 / sqrt(2 mean)`
#' (`nullb`); the simulated null gives its distribution. See
#' `proof("cconc")`.
#'
#' Given the place of each event in order (`x` as labels), the counts are
#' per place, and `grow` follows the number of distinct places as events
#' accumulate against preferential allocation, where an event goes to a
#' new place with chance `M / (M + i)` (the Chinese restaurant process):
#' `mhat` is the `M` that matches the final count, `expd` the expected
#' curve, `lo` and `hi` its bounds, and `slope` the log-log slope over the
#' second half.
#'
#' @param x The count of each unit (numbers from 0 with a positive total),
#'   or the place of each event in the order they happened (labels: text or
#'   a factor).
#' @param w `NULL`, or a multiplicity for each count (positive numbers): a
#'   banded table of people by number of placements gives each band's count
#'   once with its number of people. Counts only.
#' @param by `NULL`, or a group for each unit (each event, for places):
#'   labels. The Gini is also given within each group.
#' @param fracs The top shares to report: shares from 0 to 1 of the units.
#' @param pop `NULL`, or the population of each unit (positive numbers;
#'   counts only): the urban-scaling exponent of the counts on the
#'   population is fitted, `log x = log pref + beta log pop`, over the units
#'   with events.
#' @param nsim The number of Poisson samples at the observed mean for the
#'   null distribution of the Gini (0 for none): a whole number. It needs a
#'   whole number of units.
#' @param seed `NULL` (the draws follow [set.seed()]) or a whole number
#'   from 0.
#' @return An object of class `cconc`, a list: `stat`, a named vector
#'   (`units`, `total`, `mean`, `var`, `disp`, `ivar`, `isd`, `zero`,
#'   `znull`, `zgap`, `gini`, `gpos`, `gnull`, `vpop`, `cvb`, `nullb`);
#'   `top` (`frac`, `units`, `share`); `lorenz` (`pop`, `share`, from the
#'   origin); `null`, `NULL` or a named vector (`nsim`, `mean`, `q95`, `p`
#'   the share of simulated Gini coefficients at least the observed one
#'   with the observed counted in, and `excess` the observed minus the
#'   simulated mean); `by`, `NULL` or a data frame (`group`, `units`,
#'   `total`, `gini`, `NA` for a group without events); `scale`, `NULL` or a
#'   named vector (`beta`, `pref`, `se`, `r2`, `used`, the first four `NA`
#'   with fewer than three units with events); `grow`, `NULL` or a list
#'   (`mhat`, `slope`, `curve` with `i`, `dist`, `expd`, `lo`, `hi`); and
#'   `places`, the labels of the places (`NULL` for counts).
#' @references
#'   Weisburd, D. (2015). The law of crime concentration and the
#'   criminology of place. *Criminology*, 53(2), 133-157.
#'   \doi{10.1111/1745-9125.12070}
#'
#'   Bernasco, W. and Steenbeek, W. (2017). More places than crimes:
#'   implications for evaluating the law of crime concentration at place.
#'   *Journal of Quantitative Criminology*, 33(3), 451-467.
#'   \doi{10.1007/s10940-016-9324-7}
#'
#'   Bettencourt, L. M. A., Lobo, J., Helbing, D., Kuhnert, C. and West,
#'   G. B. (2007). Growth, innovation, scaling, and the pace of life in
#'   cities. *Proceedings of the National Academy of Sciences*, 104(17),
#'   7301-7306. \doi{10.1073/pnas.0610172104}
#' @examples
#' set.seed(2)
#' n <- rnbinom(200, mu = 2, size = 0.4)
#' r <- cconc(n, nsim = 200, seed = 1)
#' r$stat[c("gini", "zero", "gpos", "gnull")]
#' r$top
#' # people by number of placements, from a banded table
#' cconc(c(1, 2, 4, 8), w = c(500, 120, 40, 9), nsim = 0)$stat["gini"]
#' # the place of each incident, in order
#' cconc(sample(letters[1:8], 60, TRUE, prob = 8:1), nsim = 0)$grow$mhat
#' @export
cconc <- function(
  x,
  w = NULL,
  by = NULL,
  fracs = c(0.01, 0.05, 0.1, 0.25, 0.5),
  pop = NULL,
  nsim = 2000,
  seed = NULL
) {
  r <- .Call(C_clab_r_cconc, x, w, by, fracs, pop, nsim, seed)
  structure(r, class = "cconc")
}

#' Tail index of a heavy-tailed count
#'
#' The exponent `alpha` of a power-law tail, `P(X = k)` proportional to
#' `k^-alpha` above a threshold `xmin`, as when a few people account for
#' most repeat placements or most steps of an offender's path are short and
#' a few very long (Clauset, Shalizi and Newman 2009). For counts the
#' maximum likelihood estimate of the discrete law is the default, its
#' normalising constant the Hurwitz zeta function `zeta(alpha, xmin)`;
#' `approx = TRUE` gives the closed form with `xmin - 1/2`, and `discrete =
#' FALSE` the continuous Hill estimator `1 + n / sum log(x / xmin)` for
#' measurements such as distances. `ks` is the Kolmogorov-Smirnov distance
#' between the tail and the fitted law: an exponent is only as good as that
#' fit. The values at or below 0 are left out.
#'
#' @param x The values: numbers from 0.
#' @param xmin `NULL` (the 50th largest value, or the smallest when there
#'   are fewer than 50) or the threshold: one positive number.
#' @param discrete `TRUE` for counts, `FALSE` for measurements.
#' @param approx `TRUE` for the closed form for counts instead of the
#'   maximum likelihood estimate.
#' @param mintl The fewest tail values to fit to: a whole number (at least
#'   2 are always needed).
#' @return An object of class `hilltl`, a list: `alpha` and its standard
#'   error `se` (the observed information for the discrete estimate,
#'   `(alpha - 1) / sqrt(n)` otherwise; both `NA` when no fit is made),
#'   `xmin`, `ntl` (the tail size), `ks`, `reliable` (`TRUE` with 50 or more
#'   tail values and an estimate inside the search interval `[1.0001, 25]`)
#'   and `method`.
#' @references
#'   Clauset, A., Shalizi, C. R. and Newman, M. E. J. (2009). Power-law
#'   distributions in empirical data. *SIAM Review*, 51(4), 661-703.
#'   \doi{10.1137/070710111}
#'
#'   Hill, B. M. (1975). A simple general approach to inference about the
#'   tail of a distribution. *The Annals of Statistics*, 3(5), 1163-1174.
#'   \doi{10.1214/aos/1176343247}
#' @examples
#' set.seed(4)
#' k <- floor(runif(400)^(-1 / 1.5))
#' hilltl(k, xmin = 1)
#' @export
hilltl <- function(x, xmin = NULL, discrete = TRUE, approx = FALSE, mintl = 3) {
  r <- .Call(C_clab_r_hilltl, x, xmin, discrete, approx, mintl)
  structure(r, class = "hilltl")
}

#' @export
print.cconc <- function(x, ...) {
  s <- x[["stat"]]
  cat(
    "cconc: ",
    format(s[["units"]]),
    " units, total ",
    format(s[["total"]]),
    "\n",
    sprintf(
      "Gini %.4f (zero share %.4f, Gini of the rest %.4f; Poisson null %.4f)\n",
      s[["gini"]],
      s[["zero"]],
      s[["gpos"]],
      s[["gnull"]]
    ),
    sep = ""
  )
  if (!is.null(x[["null"]])) {
    cat(sprintf(
      "simulated null: mean %.4f, 95%% point %.4f, p = %.4g\n",
      x[["null"]][["mean"]],
      x[["null"]][["q95"]],
      x[["null"]][["p"]]
    ))
  }
  print(x[["top"]], row.names = FALSE)
  invisible(x)
}

#' @export
print.hilltl <- function(x, ...) {
  cat(
    "hilltl: ",
    x[["method"]],
    "\n",
    sprintf(
      "alpha %.4g (se %.3g), xmin %s, %d tail values, KS %.3g%s\n",
      x[["alpha"]],
      x[["se"]],
      format(x[["xmin"]]),
      x[["ntl"]],
      x[["ks"]],
      if (x[["reliable"]]) "" else "; not reliable"
    ),
    sep = ""
  )
  invisible(x)
}
