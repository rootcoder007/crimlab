#' Severity-weighted crime index
#'
#' Police-reported crime weighted by seriousness, as the Crime Severity
#' Index of Statistics Canada is built (Wallace et al. 2009): each offence's
#' count is multiplied by its weight (there, the incarceration rate of the
#' offence times its mean prison sentence in days), the products are summed
#' for each unit (a city, a neighbourhood) and period, divided by the
#' population and expressed against a base cell as 100.
#'
#' Give one row per unit, period and offence (rows that repeat a unit,
#' period and offence add up); an offence with no row in a cell counts 0
#' there. Beside the weighted rate the plain count rate shows whether the
#' mix of offences moved, not only their number.
#'
#' @param x The counts: numbers from 0, one per row.
#' @param off The offence of each row: labels.
#' @param w The weight of each offence: numbers from 0, named by offence
#'   (names of offences that do not occur are ignored) or in the order of
#'   the offences (a factor's levels, or the sorted labels).
#' @param unit `NULL` or the unit of each row: labels.
#' @param time `NULL` or the period of each row: labels.
#' @param pop `NULL` or the population of each row's unit and period, the
#'   same on all its rows: positive numbers.
#' @param per The rates are per `per` people: one positive number.
#' @param base `NULL`, or the row of the result whose index is 100.
#' @param sent `NULL`, or the mean sentence of each offence in the same form
#'   as `w`; the weights are then `w * sent` (an incarceration rate times a
#'   mean sentence length).
#' @return An object of class `crsev`, a list: `tab`, a data frame with one
#'   row per unit and period (in the order of the units, then the periods):
#'   `unit` and `time` when given, `count` (the counts summed), `sev` (the
#'   weighted sum), with `pop` also `pop`, `rate` and `srate` (the count and
#'   the weighted sum per `per` people), and with `base` also `idx` (`srate`,
#'   or `sev` without `pop`, as 100 times its ratio to the base row's); `w`,
#'   the weights used (`off`, `w`); `per` and `base`.
#' @references
#'   Wallace, M., Turner, J., Matarazzo, A. and Babyak, C. (2009).
#'   *Measuring Crime in Canada: Introducing the Crime Severity Index and
#'   Improvements to the Uniform Crime Reporting Survey*. Statistics
#'   Canada, catalogue 85-004-X.
#' @examples
#' d <- data.frame(city = rep(c("A", "B"), each = 2),
#'                 off = rep(c("theft", "robbery"), 2),
#'                 n = c(10, 2, 8, 3), pop = 1000)
#' crsev(d$n, d$off, c(theft = 1, robbery = 50), unit = d$city,
#'       pop = d$pop, base = 1)
#' @export
crsev <- function(
  x,
  off,
  w,
  unit = NULL,
  time = NULL,
  pop = NULL,
  per = 1e5,
  base = NULL,
  sent = NULL
) {
  r <- .Call(C_clab_r_crsev, x, off, w, unit, time, pop, per, base, sent)
  structure(r, class = "crsev")
}

#' Composite index of standardised counts
#'
#' Ranks units (neighbourhoods, divisions) on several kinds of crime at
#' once: each category's counts are standardised across the units (minus
#' the mean, over the standard deviation) and the standardised values are
#' summed with weights. A category whose counts do not vary contributes 0
#' to every unit. A unit with no row for a category counts 0 there, and
#' its standardised value is computed with that 0.
#'
#' @param x The counts: numbers from 0, one per row.
#' @param unit The unit of each row: labels, at least two units.
#' @param cat The category of each row: labels.
#' @param w `NULL` for equal weights, or the weight of each category: numbers
#'   from 0, named by category or in the order of the categories.
#' @return An object of class `zcomp`, a list: `tab`, a data frame with one
#'   row per unit (`unit`, `comp` the composite, `rank` with 1 the highest
#'   and ties sharing the lowest rank); `cat`, one row per category (`cat`,
#'   `mean`, `sd` and `w`); `z`, the matrix of standardised counts, units by
#'   categories.
#' @examples
#' d <- data.frame(hood = rep(c("a", "b", "c"), 2),
#'                 kind = rep(c("assault", "theft"), each = 3),
#'                 n = c(5, 9, 1, 30, 22, 41))
#' zcomp(d$n, d$hood, d$kind)
#' @export
zcomp <- function(x, unit, cat, w = NULL) {
  structure(.Call(C_clab_r_zcomp, x, unit, cat, w), class = "zcomp")
}

#' @export
print.crsev <- function(x, ...) {
  cat(
    "crsev: severity-weighted counts",
    if ("srate" %in% names(x[["tab"]])) {
      paste(" and rates per", format(x[["per"]]))
    },
    if (!is.na(x[["base"]])) paste(", index 100 at row", x[["base"]]),
    "\n",
    sep = ""
  )
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}

#' @export
print.zcomp <- function(x, ...) {
  cat("zcomp: composite of", nrow(x[["cat"]]), "standardised categories\n")
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}
