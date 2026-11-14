#' Segregation placements under the Mandela Rules
#'
#' The United Nations Standard Minimum Rules for the Treatment of
#' Prisoners (the Mandela Rules) define solitary confinement as 22 hours
#' or more a day without meaningful human contact (Rule 44), and prolonged
#' solitary confinement, which Rule 43 prohibits, as more than 15
#' consecutive days of it. `mndla()` classifies each placement and gives
#' the share that meets the rules, overall and by group (a year, a
#' region), with three denominators that answer different questions
#' (Sprott and Doob 2021):
#'
#' * `stays`: the share of placements;
#' * `people`: the share of people with at least one qualifying
#'   placement;
#' * `cumul`: the share of people whose days in segregation add up to more
#'   than `thr`.
#'
#' With only the length of each stay, a stay qualifies when it is longer
#' than `thr` days, and with `alert`, also when it carries at least `amin`
#' alerts (mental health, suicide risk). With `hours` (the average hours
#' out of cell a day) and `miss` (the percentage of days without the four
#' hours out of cell the law requires), the rule of Sprott and Doob: a
#' stay with at most 2 hours out and every day missed is solitary
#' confinement (class 1) up to `thr` days and torture (class 2) beyond.
#' `elig` restricts the qualifying stays (to stays with an alert, say, as a
#' proxy for the absence of contact) and `excl` removes those where
#' meaningful contact was recorded. Rows of an aggregate table enter with
#' their counts as `w`; the people denominators then need ids.
#'
#' @param days The length of each stay in days: numbers from 0.
#' @param id `NULL`, or the person of each stay: labels.
#' @param by `NULL`, or the group of each stay: labels.
#' @param w `NULL`, or the count each row stands for: numbers from 0.
#' @param thr The days beyond which a stay is prolonged: one number from 0.
#' @param hours,miss `NULL`, or the average hours out of cell (from 0) and
#'   the percentage of days without four hours out (0 to 100) of each stay.
#' @param alert `NULL`, or alerts on each stay: one column or several (a
#'   matrix, data frame or list), each `TRUE`/`FALSE`, 0/1 or
#'   `"yes"`/`"no"`.
#' @param amin The alerts at which a stay qualifies whatever its length:
#'   a whole number from 1.
#' @param elig,excl `NULL`, or for each stay `TRUE` when it may qualify
#'   (`elig`) or when contact meeting the rules was recorded (`excl`), in
#'   the forms `alert` takes.
#' @param level The level of the exact intervals: one number between 0 and
#'   1.
#' @return An object of class `mndla`, a list: `tab`, a data frame with
#'   one row per group (then `"all"`) and denominator (`group`, `den`, `n`,
#'   `k` those that qualify, `share`, `lo` and `hi` the Clopper-Pearson
#'   interval for whole counts, and with `hours` `sol` and `tor`, the
#'   stays in each class); `cls`, the class of each stay (0, 1 or 2);
#'   `alerts`, the count of alerts on each stay (with `alert`); and
#'   `rule`, `"duration"` or `"rule44"`.
#' @references
#'   United Nations (2015). *United Nations Standard Minimum Rules for the
#'   Treatment of Prisoners (the Nelson Mandela Rules)*. General Assembly
#'   resolution 70/175.
#'   <https://docs.un.org/A/RES/70/175>
#'
#'   Sprott, J. B. and Doob, A. N. (2021). *Solitary Confinement, Torture,
#'   and Canada's Structured Intervention Units*. University of Toronto,
#'   Centre for Criminology and Sociolegal Studies.
#' @examples
#' set.seed(6)
#' d <- data.frame(id = sample(1:80, 150, TRUE), yr = sample(2021:2022, 150,
#'                 TRUE), days = rgeom(150, 0.08) + 1)
#' mndla(d$days, d$id, d$yr)$tab
#' # the rule of Rules 43 and 44 with hours out of cell
#' mndla(c(5, 20, 30), hours = c(1, 1.5, 3), miss = c(100, 100, 40))$tab
#' @export
mndla <- function(
  days,
  id = NULL,
  by = NULL,
  w = NULL,
  thr = 15,
  hours = NULL,
  miss = NULL,
  alert = NULL,
  amin = 2,
  elig = NULL,
  excl = NULL,
  level = 0.95
) {
  r <- .Call(
    C_clab_r_mndla,
    days,
    id,
    by,
    w,
    thr,
    hours,
    miss,
    alert,
    amin,
    elig,
    excl,
    level
  )
  structure(r, class = "mndla")
}

#' @export
print.mndla <- function(x, ...) {
  cat(
    "mndla: placements under the Mandela Rules (",
    if (x[["rule"]] == "duration") "duration" else "Rules 43 and 44",
    ")\n",
    sep = ""
  )
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}
