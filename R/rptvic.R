#' Repeat victimisation: how often the same targets are victimised again
#'
#' Asks how much crime falls on targets (people, households, addresses,
#' businesses) already victimised, and how soon a repeat follows. From one
#' record per incident with the target's identifier (and its date):
#'
#' * incidence, prevalence and concentration (Trickett, Osborn, Seymour
#'   and Pease 1992): incidents per target at risk, victims per target at
#'   risk, and incidents per victim; and the repeats, the incidents beyond
#'   each victim's first, with their share of all incidents;
#' * the victims by number of victimisations, against the numbers a Poisson
#'   count with the same mean over the `pop` targets at risk would give:
#'   more victims with two or more than the Poisson expects (`ratio` above
#'   1) is victimisation concentrated beyond chance, the classic finding
#'   for burglary;
#' * the time between a target's successive incidents: their number, mean
#'   and median, and with `brk` the share of repeats falling in each
#'   interval since the previous incident, which for burglary is highest
#'   soon after (Johnson 2008).
#'
#' A count of victims misses the targets never victimised: `pop` supplies
#' them. The Poisson comparison assumes every target has the same risk, so
#' it measures how far risk varies, and does not by itself tell a target's
#' lasting attractiveness ("flag") from a first incident raising the risk of
#' the next ("boost").
#'
#' @param id The target of each incident: labels.
#' @param t `NULL` or the time of each incident: numbers (days, say).
#' @param pop `NULL` or the number of targets at risk: one number, at least
#'   the number of victims.
#' @param brk `NULL` or the edges of the intervals of time since the
#'   previous incident: increasing numbers, two at least (needs `t`).
#' @return An object of class `rptvic`, a list: `tab`, one row per number
#'   of victimisations `k` (from 0 with `pop`): `victims`, `incidents` and,
#'   with `pop`, `expected`, the Poisson's victims; `stat`, a named vector:
#'   `incidents`, `victims`, `repeats`, `rshare`, `conc`, `incid` and `prev`
#'   (`NA` without `pop`), `rvict` (victims with two or more), `erep` (the
#'   Poisson's number of them) and `ratio`; `gap`, `NULL` without `t`,
#'   otherwise `n`, `mean` and `median` of the times between successive
#'   incidents; and `lag`, `NULL` without `brk`, otherwise one row per
#'   interval (`lo`, `hi`, `n`, `share`, `cum`).
#' @references
#'   Trickett, A., Osborn, D. R., Seymour, J. and Pease, K. (1992). What is
#'   different about high crime areas? *The British Journal of
#'   Criminology*, 32(1), 81-89. \doi{10.1093/oxfordjournals.bjc.a048181}
#'
#'   Johnson, S. D. (2008). Repeat burglary victimisation: a tale of two
#'   theories. *Journal of Experimental Criminology*, 4(3), 215-240.
#'   \doi{10.1007/s11292-008-9055-3}
#'
#'   Farrell, G. and Pease, K. (1993). *Once Bitten, Twice Bitten: Repeat
#'   Victimisation and its Implications for Crime Prevention*. Crime
#'   Prevention Unit Series Paper 46. Home Office, London.
#' @examples
#' set.seed(8)
#' hh <- sample(1:400, 300, TRUE, prob = rgamma(400, 0.5))
#' day <- sort(sample(1:365, 300, TRUE))
#' rptvic(hh, t = day, pop = 1000, brk = c(0, 7, 30, 90, 365))
#' @export
rptvic <- function(id, t = NULL, pop = NULL, brk = NULL) {
  .Call(C_clab_r_rptvic, id, t, pop, brk)
}

#' @export
print.rptvic <- function(x, ...) {
  s <- x[["stat"]]
  cat(
    "rptvic: ",
    s[["incidents"]],
    " incidents, ",
    s[["victims"]],
    " victims, ",
    s[["repeats"]],
    " repeats (",
    format(signif(100 * s[["rshare"]], 3)),
    "%)\n",
    sep = ""
  )
  print(x[["tab"]], row.names = FALSE)
  if (!is.null(x[["lag"]])) {
    print(x[["lag"]], row.names = FALSE)
  }
  invisible(x)
}
