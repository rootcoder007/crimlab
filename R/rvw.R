#' Was the review held, and in time?
#'
#' A stay in a restrictive unit is to be reviewed by a set day. `rvwtim()`
#' takes the length of each stay and the day of its first review, and gives
#' how many stays of each length were never reviewed, how many of those past
#' the day a review was due had none, and how long the reviews that were
#' held took. Sprott, Doob and Iftene (2021, Tables 2 and 15) asked this of
#' the reviews by Independent External Decision Makers of stays in Canada's
#' Structured Intervention Units; the question is the same for any review a
#' rule sets a day for (a bail review, a segregation review, a hearing).
#'
#' A stay still open is entered with its length to the day the data end.
#' The shares carry exact (Clopper-Pearson) intervals.
#'
#' @param days The length of each stay in days: numbers from 0.
#' @param first The day of the stay on which it was first reviewed, `NA` for
#'   a stay never reviewed: numbers from 0.
#' @param due The day by which a review is due: one number above 0.
#' @param by `NULL`, or the group of each stay: labels.
#' @param cuts `NULL` (two bands of length, up to `due` and over it) or the
#'   lengths that part the bands: increasing numbers above 0.
#' @param level The level of the intervals: one number between 0 and 1.
#' @return A list: `cov`, a data frame with one row per group (then
#'   `"all"`) and band of length (`group`, `band`, `n`, `rev` the stays
#'   reviewed, `never`, `share` the share never reviewed, `lo`, `hi`);
#'   `sum`, one row per group (`group`, `n`, `rev`, `past` the stays longer
#'   than `due`, `pastnever` those of them never reviewed, `share`, `lo`,
#'   `hi`, `late` the stays first reviewed after `due`, `med` the median
#'   day of the first review); `due`; and `level`.
#' @references
#'   Sprott, J. B., Doob, A. N. and Iftene, A. (2021). *Do Independent
#'   External Decision Makers Ensure that "An Inmate's Confinement in a
#'   Structured Intervention Unit Is to End as Soon as Possible"?*
#'   Schulich School of Law, Dalhousie University, 10 May 2021.
#' @examples
#' days <- c(12, 40, 75, 90, 130, 150, 64, 20)
#' first <- c(NA, 31, 48, NA, 61, NA, 33, NA)
#' rvwtim(days, first, due = 60, cuts = c(60, 120))
#' @export
rvwtim <- function(days, first, due, by = NULL, cuts = NULL, level = 0.95) {
  fn <- "rvwtim"
  n <- .sdnum(days, "days", fn)
  if (length(first) != n || !(is.numeric(first) || all(is.na(first))) ||
        any(first < 0, na.rm = TRUE)) {
    stop("rvwtim: `first` must be numbers from 0 (NA: never reviewed), ",
         "one per stay.", call. = FALSE)
  }
  if (!.pnum1(due, .Machine$double.xmin)) {
    stop("rvwtim: `due` must be one number above 0.", call. = FALSE)
  }
  if (is.null(cuts)) cuts <- due
  .sdcuts(cuts, fn, Inf)
  g <- .sdgrp(by, n, fn)
  a <- 1 - .sdlev(level, fn)
  first <- as.numeric(first)
  rev <- !is.na(first)
  lab <- .sdband(cuts)
  bnd <- 1L + rowSums(outer(days, cuts, ">"))
  cov <- do.call(rbind, lapply(g$lev, function(k) {
    s <- g$all(k)
    nn <- tabulate(bnd[s], length(lab))
    rv <- tabulate(bnd[s & rev], length(lab))
    ci <- .sdcp(nn - rv, nn, a)
    data.frame(group = k, band = lab, n = nn, rev = rv, never = nn - rv,
               share = ci$share, lo = ci$lo, hi = ci$hi,
               stringsAsFactors = FALSE)
  }))
  sm <- do.call(rbind, lapply(g$lev, function(k) {
    s <- g$all(k)
    past <- s & days > due
    ci <- .sdcp(sum(past & !rev), sum(past), a)
    data.frame(group = k, n = sum(s), rev = sum(s & rev), past = sum(past),
               pastnever = sum(past & !rev), share = ci$share, lo = ci$lo,
               hi = ci$hi, late = sum(s & rev & first > due, na.rm = TRUE),
               med = if (any(s & rev)) stats::median(first[s & rev]) else
                 NA_real_, stringsAsFactors = FALSE)
  }))
  list(cov = cov, sum = sm, due = due, level = level)
}

#' Days on which an entitlement was missed
#'
#' The law may give a person held in a restrictive unit something every
#' day: hours out of the cell, time with other people. `entcmp()` takes,
#' for each stay, the percentage of its days on which that was missed, and
#' gives the stays by how often it was missed, the share that missed it
#' every day, and the share of all days missed. With the hours received on
#' the days it was missed, it shows how far short those days fell; with the
#' percentage of days the person refused, how many stays missed more days
#' than refusal accounts for. These are the tables of Sprott and Doob
#' (2021, Tables 6, 7, 14, 17 and 21) on the four hours out of cell and the
#' two hours of human contact in Canada's Structured Intervention Units.
#'
#' @param miss The percentage of the days of each stay on which the
#'   entitlement was missed: numbers from 0 to 100.
#' @param days `NULL`, or the length of each stay in days (numbers from 0),
#'   for the share of all days missed.
#' @param hours `NULL`, or the hours received on average on the days the
#'   entitlement was missed (numbers from 0; `NA` where none was missed).
#' @param ref `NULL`, or the percentage of the days of each stay on which
#'   the person refused: numbers from 0 to 100.
#' @param cuts The percentages that part the bands between never and
#'   always: increasing numbers above 0 and below 100.
#' @inheritParams rvwtim
#' @return A list: `tab`, a data frame with one row per group (then
#'   `"all"`) and band (`group`, `band`, `n`, `share`); `sum`, one row per
#'   group (`group`, `n`, `none` the stays that never missed, `all` those
#'   that missed every day, `share` their share, `lo`, `hi`, `half` the
#'   stays that missed on half their days or more, `dayshr` the share of
#'   all days missed, `NA` without `days`); `hrs` (with `hours`), the stays
#'   that missed by the hours they received (`band`, `n`, `share`); `ref`
#'   (with `ref`), one row per group (`group`, `n`, `most` the stays that
#'   refused on more than half their days, `net` those that missed on more
#'   days than they refused); and `level`.
#' @references
#'   Sprott, J. B. and Doob, A. N. (2021). *Solitary Confinement, Torture,
#'   and Canada's Structured Intervention Units*. University of Toronto,
#'   Centre for Criminology and Sociolegal Studies, 23 February 2021.
#' @examples
#' miss <- c(100, 100, 80, 45, 10, 0, 100, 60)
#' entcmp(miss, days = c(20, 35, 16, 40, 18, 5, 62, 30),
#'        hours = c(0.4, 1.5, 2.2, 3, 3.5, NA, 0.8, 1),
#'        ref = c(10, 100, 5, 0, 0, 0, 30, 60))
#' @export
entcmp <- function(miss, days = NULL, by = NULL, hours = NULL, ref = NULL,
                   cuts = c(20, 50), level = 0.95) {
  fn <- "entcmp"
  n <- .sdpct(miss, "miss", fn)
  if (!is.null(days) && .sdnum(days, "days", fn) != n) .sdlen("days", fn)
  .sdcuts(cuts, fn, 100)
  g <- .sdgrp(by, n, fn)
  a <- 1 - .sdlev(level, fn)
  lab <- c("0", .sdband(cuts, 100), "100")
  bnd <- 2L + rowSums(outer(miss, cuts, ">"))
  bnd[miss == 0] <- 1L
  bnd[miss == 100] <- length(lab)
  tab <- do.call(rbind, lapply(g$lev, function(k) {
    s <- g$all(k)
    nn <- tabulate(bnd[s], length(lab))
    data.frame(group = k, band = lab, n = nn, share = nn / sum(s),
               stringsAsFactors = FALSE)
  }))
  sm <- do.call(rbind, lapply(g$lev, function(k) {
    s <- g$all(k)
    ci <- .sdcp(sum(s & miss == 100), sum(s), a)
    data.frame(group = k, n = sum(s), none = sum(s & miss == 0),
               all = sum(s & miss == 100), share = ci$share, lo = ci$lo,
               hi = ci$hi, half = sum(s & miss >= 50),
               dayshr = if (is.null(days)) NA_real_ else
                 sum(miss[s] * days[s]) / (100 * sum(days[s])),
               stringsAsFactors = FALSE)
  }))
  out <- list(tab = tab, sum = sm)
  if (!is.null(hours)) out$hrs <- .enthrs(hours, miss, n)
  if (!is.null(ref)) out$ref <- .entref(ref, miss, n, g)
  out$level <- level
  out
}

#' Do decision makers decide alike?
#'
#' When cases are shared among several decision makers (reviewers, judges,
#' panels), the share of their decisions that go one way should differ by
#' no more than chance if the cases are alike. `dmvar()` gives each
#' decision maker's rate with its exact interval and the share of all such
#' decisions that are theirs, the chi-square test that the rates are equal,
#' the same test without the decision makers who decided few cases, and the
#' test with each decision maker left out in turn, which shows whether one
#' of them makes the difference. Sprott, Doob and Iftene (2021, Tables 10
#' and 11) found the share of reviews that kept a person in a Structured
#' Intervention Unit to run from 37.5% to over 75% across twelve decision
#' makers, with and without the one who had decided fewest.
#'
#' The test takes the cases as shared without regard to who decides; where
#' they are not, a difference may be one of the cases and not of the
#' decision makers ([jdgiv()] and [jmono()] treat that question).
#'
#' @param dec The decision in each case: labels.
#' @param dm The decision maker of each case: labels.
#' @param out The decision whose rate is compared: one of the labels of
#'   `dec`.
#' @param nmin The cases below which a decision maker is left out of the
#'   second test: a whole number from 1.
#' @inheritParams rvwtim
#' @return A list: `tab`, a data frame with one row per decision maker
#'   (`dm`, `n`, `k` the decisions that were `out`, `rate`, `lo`, `hi`,
#'   `share` of all such decisions); `test`, two rows (`set`: `"all"` and
#'   `"n >= nmin"`, with `m` the decision makers, `stat`, `df`, `p` and
#'   `low` the cells with an expected count under 5); `drop`, one row per
#'   decision maker left out (`dm`, `stat`, `df`, `p`); `out`; and `level`.
#' @references
#'   Sprott, J. B., Doob, A. N. and Iftene, A. (2021). *Do Independent
#'   External Decision Makers Ensure that "An Inmate's Confinement in a
#'   Structured Intervention Unit Is to End as Soon as Possible"?*
#'   Schulich School of Law, Dalhousie University, 10 May 2021.
#' @examples
#' dm <- rep(c("a", "b", "c", "d"), c(30, 28, 25, 4))
#' dec <- rep(c("remain", "remove", "remain", "remove", "remain", "remove",
#'              "remain"), c(24, 6, 10, 18, 15, 10, 4))
#' dmvar(dec, dm, out = "remain")
#' @export
dmvar <- function(dec, dm, out, nmin = 5, level = 0.95) {
  .dmck(dec, dm, out, nmin)
  dec <- as.character(dec)
  dm <- as.character(dm)
  lev <- unique(dm)
  a <- 1 - .sdlev(level, "dmvar")
  j <- match(dm, lev)
  nn <- tabulate(j, length(lev))
  kk <- tabulate(j[dec == as.character(out)], length(lev))
  ci <- .sdcp(kk, nn, a)
  tab <- data.frame(dm = lev, n = nn, k = kk, rate = ci$share, lo = ci$lo,
                    hi = ci$hi, share = kk / sum(kk),
                    stringsAsFactors = FALSE)
  big <- nn >= nmin
  test <- cbind(set = c("all", "n >= nmin"),
                rbind(.dmchi(nn, kk), .dmchi(nn[big], kk[big])),
                stringsAsFactors = FALSE)
  drop <- do.call(rbind, lapply(seq_along(lev), function(i) {
    cbind(dm = lev[i], .dmchi(nn[-i], kk[-i])[c("stat", "df", "p")],
          stringsAsFactors = FALSE)
  }))
  list(tab = tab, test = test, drop = drop, out = as.character(out),
       level = level)
}

# ---- checks and parts shared by the three --------------------------------

.sdnum <- function(x, nm, fn) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(x < 0)) {
    stop(fn, ": `", nm, "` must be numbers from 0 with no missing value.",
         call. = FALSE)
  }
  length(x)
}

.sdpct <- function(x, nm, fn) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(x < 0 | x > 100)) {
    stop(fn, ": `", nm, "` must be percentages from 0 to 100 with no ",
         "missing value.", call. = FALSE)
  }
  length(x)
}

.sdlen <- function(nm, fn) {
  stop(fn, ": `", nm, "` must have one value per stay.", call. = FALSE)
}

.sdcuts <- function(cuts, fn, top) {
  bad <- !is.numeric(cuts) || !length(cuts) || anyNA(cuts)
  if (!bad) {
    bad <- any(cuts <= 0 | cuts >= top) || is.unsorted(cuts, strictly = TRUE)
  }
  if (bad) {
    stop(fn, ": `cuts` must be increasing numbers above 0",
         c("", paste(" and below", top))[is.finite(top) + 1L], ".",
         call. = FALSE)
  }
  invisible(NULL)
}

.sdlev <- function(level, fn) {
  if (!.pnum1(level) || level <= 0 || level >= 1) {
    stop(fn, ": `level` must be one number between 0 and 1.", call. = FALSE)
  }
  level
}

# the groups in the order they appear, then "all"
.sdgrp <- function(by, n, fn) {
  if (is.null(by)) {
    return(list(lev = "all", all = function(k) rep(TRUE, n)))
  }
  if (length(by) != n || anyNA(by) || !is.atomic(by)) {
    stop(fn, ": `by` must have one label per stay and no missing value.",
         call. = FALSE)
  }
  by <- as.character(by)
  list(lev = c(unique(by), "all"),
       all = function(k) if (k == "all") rep(TRUE, n) else by == k)
}

# the names of the bands the cuts part: up to the first, between, over the
# last (or, with `top`, over the last and under `top`)
.sdband <- function(cuts, top = NULL) {
  f <- sprintf("%.15g", cuts)
  k <- length(f)
  c(if (is.null(top)) paste("up to", f[1L]) else paste("over 0 to", f[1L]),
    if (k > 1L) paste("over", f[-k], "to", f[-1L]),
    if (is.null(top)) paste("over", f[k]) else
      paste("over", f[k], "to under", top))
}

# k of n: the share and its Clopper-Pearson limits (NA of no cases)
.sdcp <- function(k, n, a) {
  lo <- suppressWarnings(ifelse(k > 0, stats::qbeta(a / 2, k, n - k + 1), 0))
  hi <- suppressWarnings(ifelse(k < n, stats::qbeta(1 - a / 2, k + 1, n - k),
                                1))
  z <- n == 0
  list(share = ifelse(z, NA_real_, k / n), lo = ifelse(z, NA_real_, lo),
       hi = ifelse(z, NA_real_, hi))
}

# the stays that missed, by the hours they received on the days they missed
.enthrs <- function(hours, miss, n) {
  if (length(hours) != n || !(is.numeric(hours) || all(is.na(hours))) ||
        any(hours < 0, na.rm = TRUE)) {
    stop("entcmp: `hours` must be numbers from 0 (NA: none missed), one ",
         "per stay.", call. = FALSE)
  }
  hc <- c(0.5, 1, 2, 3)
  s <- miss > 0 & !is.na(hours)
  nn <- tabulate(1L + rowSums(outer(as.numeric(hours)[s], hc, ">")), 5L)
  data.frame(band = .sdband(hc), n = nn, share = nn / sum(s),
             stringsAsFactors = FALSE)
}

# refusals: the stays that refused most days, and those that missed on more
# days than they refused
.entref <- function(ref, miss, n, g) {
  if (.sdpct(ref, "ref", "entcmp") != n) .sdlen("ref", "entcmp")
  do.call(rbind, lapply(g$lev, function(k) {
    s <- g$all(k)
    data.frame(group = k, n = sum(s), most = sum(s & ref > 50),
               net = sum(s & miss > ref), stringsAsFactors = FALSE)
  }))
}

# the checks of dmvar()
.dmck <- function(dec, dm, out, nmin) {
  n <- length(dec)
  if (!n || anyNA(dec) || !is.atomic(dec)) {
    stop("dmvar: `dec` must be labels with no missing value.", call. = FALSE)
  }
  if (length(dm) != n || anyNA(dm) || !is.atomic(dm)) {
    stop("dmvar: `dm` must have one label per case and no missing value.",
         call. = FALSE)
  }
  if (length(out) != 1L || is.na(out) ||
        !as.character(out) %in% as.character(dec)) {
    stop("dmvar: `out` must be one of the decisions in `dec`.", call. = FALSE)
  }
  if (length(unique(as.character(dm))) < 2L) {
    stop("dmvar: `dm` must name two or more decision makers.", call. = FALSE)
  }
  if (!.pnum1(nmin, 1, Inf, TRUE)) {
    stop("dmvar: `nmin` must be one whole number from 1.", call. = FALSE)
  }
  invisible(NULL)
}

# Pearson's chi-square that k of n is one rate over the decision makers
.dmchi <- function(n, k) {
  m <- length(n)
  tot <- sum(n)
  p <- sum(k) / tot
  if (m < 2L || p == 0 || p == 1) {
    return(data.frame(m = m, stat = NA_real_, df = NA_real_, p = NA_real_,
                      low = NA_real_))
  }
  e <- c(n * p, n * (1 - p))
  st <- sum((c(k, n - k) - e)^2 / e)
  data.frame(m = m, stat = st, df = as.numeric(m - 1),
             p = stats::pchisq(st, m - 1, lower.tail = FALSE),
             low = as.numeric(sum(e < 5)))
}
