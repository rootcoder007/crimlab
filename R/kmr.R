#' Kaplan-Meier survival and recurrence
#'
#' Estimates the share still without an event at each time since the
#' start (release, placement, a case opening, the last incident at a
#' place) by Kaplan-Meier, by group, with Greenwood's standard error and
#' the log-scale interval of `survival::survfit()` (its default), each
#' group's median (where the curve sits at exactly 1/2, halfway to its
#' next fall, as `survfit()` does), and the log-rank test between groups
#' (`survival::survdiff()`). Each product step is a factor in \eqn{[0, 1]},
#' so the curve starts at most 1 and never rises (`proof("kmrec")`).
#'
#' Durations are given directly (`time`, with `event` 0 where a duration
#' was cut short), or, with `unit`, `time` holds the times of events and
#' the durations are the gaps between one unit's consecutive events (time
#' to a repeat incident at a place, a repeat admission of a person); with
#' `end` each unit's last gap, from its last event to `end`, enters as
#' censored, so units with long quiet spells count.
#'
#' @param time Durations from 0, or with `unit` the times of events:
#'   finite numbers from 0.
#' @param event `NULL` (every duration ends in an event) or 1 for an event
#'   and 0 for a censored duration (logical or numbers); not read with
#'   `unit`.
#' @param grp `NULL`, or the group of each time: labels. With `unit` it
#'   must be the same for every time of a unit.
#' @param unit `NULL`, or the unit of each event time: labels.
#' @param end With `unit`, `NULL` or the end of observation: one number at
#'   or after every time.
#' @param level The confidence level: one number strictly between 0 and 1.
#' @return An object of class `kmrec`, a list: `km`, one row per group and
#'   distinct time (`grp`, `time`, `nrisk`, `nevent`, `ncens`, `surv`,
#'   `se`, `lo`, `hi`); `med`, one row per group (`grp`, `n`, `events`,
#'   `median`, `NA` when the curve never reaches 1/2); `lr`, `NULL` for one
#'   group, otherwise the log-rank `stat`, `df` and `p`; `n`, the durations;
#'   and `level`.
#' @references
#'   Kaplan, E. L. and Meier, P. (1958). Nonparametric estimation from
#'   incomplete observations. *Journal of the American Statistical
#'   Association*, 53(282), 457-481. \doi{10.1080/01621459.1958.10501452}
#'
#'   Maltz, M. D. (1984). *Recidivism*. Academic Press.
#' @examples
#' kmrec(c(3, 8, 8, 12, 20, 24), event = c(1, 1, 0, 1, 0, 1),
#'       grp = c("a", "a", "a", "b", "b", "b"))
#' # repeat incidents at three places, watched to day 100
#' kmrec(c(2, 9, 30, 5, 70, 11, 12, 13),
#'       unit = rep(c("p", "q", "r"), c(3, 2, 3)), end = 100)$med
#' @export
kmrec <- function(
  time,
  event = NULL,
  grp = NULL,
  unit = NULL,
  end = NULL,
  level = 0.95
) {
  .Call(C_clab_r_kmrec, time, event, grp, unit, end, level)
}

#' @export
print.kmrec <- function(x, ...) {
  cat("kmrec: ", x[["n"]], " durations\n", sep = "")
  print(x[["med"]], row.names = FALSE)
  l <- x[["lr"]]
  if (!is.null(l)) {
    cat("log-rank chi-square ", format(l[["stat"]]), " on ", l[["df"]],
        " df, p ", format(l[["p"]]), "\n", sep = "")
  }
  invisible(x)
}
