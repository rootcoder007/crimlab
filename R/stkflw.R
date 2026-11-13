#' Stock and flow of a population held over time
#'
#' Person-days divided by the length of a period give the average daily
#' population (a stock: how many are held at once); divided by the people
#' admitted or released, the average length of stay (how long each
#' stays). The two share a numerator and `days = adp * t = people * alos`
#' (Little's law in the words of corrections planning; Lakner 1976), so
#' any two of people, days, average daily population and average length
#' of stay give the rest, except days with the average daily population,
#' which carry the same information. When the length of stay moves, the
#' stock and the flow can move in opposite directions; `pshare` is the
#' share of the change in days (on the log scale) that the change in
#' people accounts for, the rest coming from the length of stay. See
#' `proof("stkflw")`.
#'
#' Given the releases of each period and a number of days by which every
#' release is brought forward, `cut` is the stock that removes in a steady
#' state, `rel * early / t` (releasing everyone six months early removes
#' six months of releases at once), and `cutpct` that share of the average
#' daily population.
#'
#' @param people,days,adp,alos Two of: the people admitted or released in
#'   each period (positive), the person-days (from 0), the average daily
#'   population (from 0; or a list with the periodic headcounts of each
#'   period, whose means are taken), and the average length of stay
#'   (positive).
#' @param t The length of each period in days: one positive number or one
#'   per period.
#' @param prd `NULL` (1, 2, ...) or the label of each period.
#' @param from,to `NULL`, or the first and last day of each period (dates
#'   or date strings); then `t` is the days from one to the other,
#'   counting both.
#' @param expo `NULL`, or the population at risk: one positive number or
#'   one per period, for rates per `per`.
#' @param per The rates are per `per` people: one positive number.
#' @param base Changes against the `"first"` period or the `"previous"`
#'   one.
#' @param rel,early `NULL`, or the releases of each period and the days by
#'   which every release is brought forward.
#' @return An object of class `stkflw`, a list: `flow`, a data frame with
#'   one row per period (`period`, `t`, `people`, `days`, `adp`, `alos`,
#'   their percent changes `peochg`, `aloschg`, `dayschg`, `adpchg`,
#'   `pshare`, with `expo` `flowrate`, `stkrate`, `frchg`, `srchg`, and
#'   with `rel` `cut`, `cutpct`); `mean`, the means over the periods of
#'   `people`, `days`, `adp` and `alos`; `base` and `per`.
#' @references
#'   Lakner, E. (1976). *A Manual of Statistical Sampling Methods for
#'   Corrections Planners*. University of Illinois at Urbana-Champaign.
#'
#'   Little, J. D. C. (1961). A proof for the queuing formula: L = lambda
#'   W. *Operations Research*, 9(3), 383-387. \doi{10.1287/opre.9.3.383}
#' @examples
#' stkflw(people = c(5200, 4800, 4300), days = c(310000, 318000, 330000),
#'        prd = 2021:2023)
#' # admissions implied by a population and a length of stay
#' stkflw(adp = 900, alos = 61)$flow$people
#' @export
stkflw <- function(
  people = NULL,
  days = NULL,
  adp = NULL,
  alos = NULL,
  t = 365,
  prd = NULL,
  from = NULL,
  to = NULL,
  expo = NULL,
  per = 1e5,
  base = "first",
  rel = NULL,
  early = NULL
) {
  if (!is.null(from) || !is.null(to)) {
    t <- .perdays(from, to)
  }
  r <- .Call(
    C_clab_r_stkflw,
    people,
    days,
    adp,
    alos,
    t,
    prd,
    expo,
    per,
    base,
    rel,
    early
  )
  structure(r, class = "stkflw")
}

#' Lengths of stay: their distribution and its fit
#'
#' The lengths of stay of a set of placements (or of people): the mean
#' with its t interval, the median and quartiles, the longest; the share
#' of stays in each band and at or below each threshold (15 days, the
#' Mandela Rules' limit, say), with the share of each band carrying a flag
#' (no review on record); and the fit of lognormal, exponential and Pareto
#' laws by maximum likelihood, compared by AIC: a heavy tail (Pareto or
#' lognormal beating the exponential) separates a few long stays from many
#' short ones. A stay still open is a lower bound on its length; see
#' [losbnd()] for the mean it bounds.
#'
#' @param days The stays: numbers from 0.
#' @param w `NULL`, or how many stays each value stands for: whole numbers
#'   from 0 with a positive total.
#' @param breaks `NULL`, or band edges, increasing; a band is `(lo, hi]`.
#' @param thr `NULL`, or thresholds: numbers from 0.
#' @param flag `NULL`, or `TRUE`/`FALSE` for each stay.
#' @param level The level of the interval of the mean: one number between
#'   0 and 1.
#' @return An object of class `lostay`, a list: `stay`, a named vector
#'   (`n`, `days`, `mean`, `sd`, `median`, `q1`, `q3`, `max`, `se`, `lo`,
#'   `hi`); `bands` (with `breaks`: `lo`, `hi`, `n`, `share` of all stays,
#'   with `flag` `flagged` the share of the band's stays); `thr` (`thr`,
#'   `atmost`, `above`); `fit` (`model`, `aic`, `best`; the fits use the
#'   stays above 0); and `par` (`meanlog`, `sdlog`, `alpha`).
#' @references Clauset, A., Shalizi, C. R. and Newman, M. E. J. (2009).
#'   Power-law distributions in empirical data. *SIAM Review*, 51(4),
#'   661-703. \doi{10.1137/070710111}
#' @examples
#' set.seed(3)
#' s <- round(rlnorm(400, 2.3, 1.1)) + 1
#' r <- lostay(s, breaks = c(0, 5, 15, 31, 61, Inf), thr = 15)
#' r$bands
#' r$fit
#' @export
lostay <- function(
  days,
  w = NULL,
  breaks = NULL,
  thr = NULL,
  flag = NULL,
  level = 0.95
) {
  r <- .Call(C_clab_r_lostay, days, w, breaks, thr, flag, level)
  structure(r, class = "lostay")
}

#' Little's law on a court docket
#'
#' From the arrival and disposition times of cases (a disposition missing
#' or after the horizon: still pending), the time-average backlog equals
#' the filing rate times the mean time to disposition when every case has
#' been disposed (Little 1961): `resid` is that identity's residual, 0 up
#' to rounding. The backlog counts each pending case to the horizon; the
#' mean time is over the disposed cases, which understates the truth when
#' cases are pending (see [losbnd()]). Given a target backlog, `need` is
#' the mean time to disposition that holds it at the present filing rate.
#' See `proof("bklog")`.
#'
#' @param arrive The arrival time of each case: numbers from 0.
#' @param dispose The disposition time of each case: numbers, `NA` for one
#'   pending.
#' @param horizon `NULL` (the latest arrival or disposition) or the end of
#'   the window: one positive number.
#' @param target `NULL`, or a target backlog: one number from 0.
#' @return An object of class `bklog`, a list: `est`, a named vector
#'   (`cases`, `pending`, `rate` the filing rate, `wait` the mean time to
#'   disposition, `backlog`, `occup` the case-time, `resid`); `horizon`;
#'   and `need` (with `target`).
#' @references Little, J. D. C. (1961). A proof for the queuing formula: L
#'   = lambda W. *Operations Research*, 9(3), 383-387.
#'   \doi{10.1287/opre.9.3.383}
#' @examples
#' bklog(arrive = c(0, 2, 3, 7, 9), dispose = c(5, 4, 12, NA, 15),
#'       horizon = 20, target = 2)
#' @export
bklog <- function(arrive, dispose, horizon = NULL, target = NULL) {
  structure(
    .Call(C_clab_r_bklog, arrive, dispose, horizon, target),
    class = "bklog"
  )
}

#' Case flow through a court
#'
#' Cases received and disposed by period: the clearance rate (disposed
#' over received, above 1 when the caseload shrinks) with its exact
#' interval, the net, the clearance over a rolling window, and with the
#' pending caseload its change and the periods it would take to clear at
#' the recent pace. With the length of each period (`span`) the filing
#' rate gives, by Little's law, the mean time to disposition the pending
#' caseload implies, to set beside a published mean (`mtime`): their
#' `ratio` far from 1 means the published mean leaves out the long-pending
#' cases. Published duration bands (`shares` of cases in bands starting
#' at `floors`) bound the mean from below.
#'
#' With `mix`, the counts of each outcome (conviction, withdrawal,
#' acquittal) by period: each period's shares with Goodman's simultaneous
#' intervals, the chi-square test that the mix is the same in every
#' period with Cramer's V, and the Cochran-Armitage trend in each share
#' over the periods (as [stats::prop.trend.test()]) with Holm's
#' adjustment.
#'
#' @param recv,disp `NULL`, or the cases received and disposed in each
#'   period: numbers from 0.
#' @param pend `NULL`, or the cases pending at the end of each period.
#' @param prd `NULL` (1, 2, ...) or the label of each period.
#' @param win The rolling window in periods: a whole number from 1.
#' @param level The level of the intervals: one number between 0 and 1.
#' @param span `NULL`, or the days in each period (one number or one per
#'   period); needs `pend`.
#' @param mtime `NULL`, or the published mean days to disposition of each
#'   period; needs `span`.
#' @param shares,floors `NULL`, or the share of cases in each duration band
#'   (a matrix, one row per period) and each band's lower edge.
#' @param mix `NULL`, or the outcome counts: a numeric matrix of periods by
#'   outcomes (named by its dimnames).
#' @param scores `NULL` (1, 2, ...) or the period scores of the trend test.
#' @return An object of class `dktflw`, a list: `flow` (with `recv`:
#'   `period`, `recv`, `disp`, `clear`, `lo`, `hi`, `net`, `rclear`,
#'   `rlo`, `rhi`, then `pend`, `pchg`, `toclear`, or `cumnet`, and with
#'   `span` `frate`, `impl`, with `mtime` `ratio`, `ltlbk`, with `shares`
#'   `floorbd`); `mix` (a list of `shares`, `homog` and `trend`); and
#'   `level`.
#' @references Goodman, L. A. (1965). On simultaneous confidence intervals
#'   for multinomial proportions. *Technometrics*, 7(2), 247-254.
#'   \doi{10.1080/00401706.1965.10490252}
#' @examples
#' dktflw(recv = c(120, 130, 125, 140), disp = c(110, 135, 120, 150),
#'        pend = c(400, 395, 400, 390), win = 2, span = 91)
#' m <- matrix(c(50, 55, 60, 30, 28, 25, 20, 17, 15), 3,
#'             dimnames = list(2021:2023, c("guilty", "withdrawn", "other")))
#' dktflw(mix = m)$mix$trend
#' @export
dktflw <- function(
  recv = NULL,
  disp = NULL,
  pend = NULL,
  prd = NULL,
  win = 12,
  level = 0.95,
  span = NULL,
  mtime = NULL,
  shares = NULL,
  floors = NULL,
  mix = NULL,
  scores = NULL
) {
  r <- .Call(
    C_clab_r_dktflw,
    recv,
    disp,
    pend,
    prd,
    win,
    level,
    span,
    mtime,
    shares,
    floors,
    mix,
    scores
  )
  structure(r, class = "dktflw")
}

#' Actuarial life table of grouped times
#'
#' How long cases take, from the counts of cases ending (and censored) in
#' each band of time: the conditional probability of ending in each band
#' with the censored counted as half a case at risk, the survival (the
#' share still open at the end of each band) with Greenwood's standard
#' error and a plain, log or log-log interval, the hazard, and the median,
#' linear within its band.
#'
#' @param breaks The band edges, increasing (the last may be `Inf`), or
#'   band labels as [bandx()] reads them (`"0-3"`, `"Over 12"`), contiguous.
#' @param events The cases ending in each band: numbers from 0.
#' @param cens The cases censored in each band: one number or one per band.
#' @param level The level of the intervals: one number between 0 and 1.
#' @param ci `"log-log"`, `"log"` or `"plain"`.
#' @return An object of class `lftab`, a list: `tab` (`lower`, `upper`,
#'   `entering`, `atrisk`, `events`, `cens`, `q`, `surv`, `se`, `lo`, `hi`,
#'   `hazard`), `median` (`Inf` when survival stays above one half or
#'   crosses it in an open band) and `ci`.
#' @examples
#' lftab(c(0, 3, 6, 12, 24, Inf), events = c(40, 30, 20, 10, 5),
#'       cens = c(2, 2, 3, 3, 4))
#' @export
lftab <- function(breaks, events, cens = 0, level = 0.95, ci = "log-log") {
  if (is.character(breaks) || is.factor(breaks)) {
    b <- .bands(breaks, TRUE, FALSE)
    lo <- b$lower
    lo[b$openlo] <- 0
    hi <- ifelse(b$openup, Inf, b$upper)
    if (anyNA(lo) || anyNA(hi)) {
      stop(
        "lftab: these bands cannot be read: ",
        toString(b$label[is.na(lo) | is.na(hi)]),
        ".",
        call. = FALSE
      )
    }
  } else {
    if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks)) {
      stop("lftab: `breaks` must be at least two band edges.", call. = FALSE)
    }
    lo <- as.numeric(breaks[-length(breaks)])
    hi <- as.numeric(breaks[-1L])
  }
  r <- .Call(C_clab_r_lftab, lo, hi, events, cens, level, ci)
  structure(r, class = "lftab")
}

#' Transition matrix between states
#'
#' The moves of people between states (regions, institutions, security
#' levels) from their sequence of records: the counts and probabilities of
#' each move between consecutive records of the same person, the
#' stationary distribution of the chain (where the population would settle
#' if the probabilities held), the share of moves that stay in their
#' state, and the Theil index of the moves between states (0 when they
#' spread evenly over the pairs of states, higher when a few pairs carry
#' them).
#'
#' @param state The state of each record: labels.
#' @param id The person (or sequence) of each record: labels. Give
#'   `paste(id, year)` to keep the moves within a year.
#' @param time `NULL` (the order of the rows) or the time of each record:
#'   numbers; ties keep the order of the rows.
#' @return An object of class `trnmat`, a list: `count` and `prob`, the
#'   matrices of moves from (rows) and to (columns); `pi`, the stationary
#'   distribution (`NA` when a state is never left or it is not unique);
#'   and `summ` (`moves`, `stay`, `theil`).
#' @examples
#' r <- trnmat(c("N", "N", "E", "N", "W", "W", "E", "E"),
#'             c(1, 1, 1, 2, 2, 2, 3, 3), time = c(1, 2, 3, 1, 2, 3, 1, 2))
#' r$prob
#' r$pi
#' @export
trnmat <- function(state, id, time = NULL) {
  structure(.Call(C_clab_r_trnmat, state, id, time), class = "trnmat")
}

#' @export
print.stkflw <- function(x, ...) {
  cat("stkflw: stock and flow, changes against the", x[["base"]], "period\n")
  print(x[["flow"]], row.names = FALSE)
  invisible(x)
}

#' @export
print.lostay <- function(x, ...) {
  s <- x[["stay"]]
  cat(sprintf(
    "lostay: %s stays, mean %.4g (%.4g to %.4g), median %.4g\n",
    format(s[["n"]]),
    s[["mean"]],
    s[["lo"]],
    s[["hi"]],
    s[["median"]]
  ))
  if (!is.null(x[["bands"]])) {
    print(x[["bands"]], row.names = FALSE)
  }
  print(x[["fit"]], row.names = FALSE)
  invisible(x)
}

#' @export
print.bklog <- function(x, ...) {
  cat("bklog: Little's law over a horizon of", format(x[["horizon"]]), "\n")
  print(x[["est"]])
  invisible(x)
}

#' @export
print.dktflw <- function(x, ...) {
  cat("dktflw: case flow\n")
  if (!is.null(x[["flow"]])) {
    print(x[["flow"]], row.names = FALSE)
  }
  if (!is.null(x[["mix"]])) {
    print(x[["mix"]][["trend"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.lftab <- function(x, ...) {
  cat(
    "lftab: life table,",
    x[["ci"]],
    "intervals, median",
    format(x[["median"]]),
    "\n"
  )
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}

#' @export
print.trnmat <- function(x, ...) {
  s <- x[["summ"]]
  cat(sprintf(
    "trnmat: %s moves, %.3g stay in their state\n",
    format(s[["moves"]]),
    s[["stay"]]
  ))
  print(x[["prob"]])
  invisible(x)
}

# Days from `from` to `to`, counting both: one per period.
.perdays <- function(from, to) {
  rd <- function(v, nm) {
    if (
      is.null(v) ||
        is.numeric(v) ||
        !(inherits(v, c("Date", "POSIXt")) ||
          is.character(v))
    ) {
      stop(
        "stkflw: `from` and `to` must both be dates or date strings.",
        call. = FALSE
      )
    }
    d <- tryCatch(as.Date(v), error = function(e) NA)
    if (anyNA(d)) {
      stop(
        "stkflw: `",
        nm,
        "` holds something that is not a date.",
        call. = FALSE
      )
    }
    d
  }
  a <- rd(from, "from")
  b <- rd(to, "to")
  if (length(a) != length(b) || any(b < a)) {
    stop("stkflw: `to` must not precede `from`, one per period.", call. = FALSE)
  }
  as.numeric(b - a) + 1
}
