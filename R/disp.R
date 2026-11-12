#' Veil-of-darkness test
#'
#' Compares each group's share of stops made in daylight, when a driver's
#' race can be seen, with its share in darkness, when it cannot. The ratio of
#' the two odds against a reference group is the factor by which visibility
#' changes that group's chance of being stopped; 1 means no change.
#'
#' The ratio identifies the visibility effect only when the mix of drivers on
#' the road is the same in daylight and darkness. `shift` states how much each
#' group's exposure differs between the two, and `adj` is the ratio with that
#' difference taken out. See `proof("vodtst")`.
#'
#' @param day Stops by group in daylight: non-negative numbers, optionally
#'   named.
#' @param night Stops by group in darkness: positive numbers, same length.
#' @param ref The reference group, by name or position.
#' @param shift `NULL`, or one positive number per group: daylight exposure
#'   over darkness exposure.
#' @return A data frame with one row per group: `group`, `share_day`,
#'   `share_night`, `ratio`, `sign` (the sign of `share_day - share_night`)
#'   and, when `shift` is given, `adj`.
#' @references Grogger, J. and Ridgeway, G. (2006). Testing for racial
#'   profiling in traffic stops from behind a veil of darkness. *Journal of
#'   the American Statistical Association*, 101(475), 878-887.
#'   \doi{10.1198/016214506000000168}
#' @examples
#' vodtst(day = c(white = 600, black = 400),
#'        night = c(white = 700, black = 300))
#' @export
vodtst <- function(day, night, ref = 1, shift = NULL) {
  nms <- .grp(day)
  day <- .num(day, "day", lo = 0)
  n <- length(day)
  night <- .num(night, "night", len = n, lo = 0, open = TRUE)
  if (!is.null(shift)) {
    shift <- .num(shift, "shift", len = n, lo = 0, open = TRUE)
  }
  i <- .ref(ref, nms, n)
  if (day[i] <= 0) {
    stop("the reference group has no daylight stops.", call. = FALSE)
  }
  r <- .Call(C_csls_vod, day, night, i - 1L, shift)
  out <- data.frame(group = nms, share_day = r[["share_day"]],
                    share_night = r[["share_night"]], ratio = r[["ratio"]],
                    sign = r[["sign"]])
  if (!is.null(shift)) out[["adj"]] <- r[["adj"]]
  out
}

#' Risks behind an odds ratio
#'
#' An odds ratio is not a relative risk: it lies further from 1, and by more
#' the more common the outcome. Given the odds ratio, the overall rate of the
#' outcome and the share exposed, this returns the two risks and the relative
#' risk they imply. See `proof("orrisk")`.
#'
#' @param or The odds ratio: one positive number.
#' @param rate The overall rate of the outcome, strictly between 0 and 1.
#' @param share The exposed share of the population, strictly between 0 and 1.
#' @return A named numeric vector: `exposed` and `unexposed` (the two risks),
#'   `rr` (their ratio) and `or`. A risk closer to 1 than a double can hold
#'   is returned as 1.
#' @examples
#' orrisk(or = 3, rate = 0.4, share = 0.5)
#' @export
orrisk <- function(or, rate, share) {
  or <- .num(or, "or", len = 1L, lo = 0, open = TRUE)
  if (!is.finite(1 / or)) {
    stop("`or` is too small: its reciprocal is not finite.", call. = FALSE)
  }
  rate <- .one(rate, "rate", 0, 1)
  share <- .one(share, "share", 0, 1)
  r <- .Call(C_csls_orrisk, or, rate, share)
  c(exposed = r[1L], unexposed = r[2L], rr = r[1L] / r[2L], or = or)
}

#' Outcome test
#'
#' The stop rate and the hit rate of a group searched whenever the chance of
#' guilt exceeds a threshold. Hit rates compare averages, not thresholds: two
#' groups held to the same threshold can show different hit rates, and a group
#' held to a lower threshold can show the higher one (infra-marginality). See
#' `proof("outtst")`.
#'
#' @param q The chances of guilt in the group: numbers in the interval 0 to 1.
#' @param t The search threshold: one number.
#' @param w `NULL` for equal weights, or one non-negative weight per element
#'   of `q` with a positive total.
#' @return A named numeric vector: `stop` (the weighted share searched) and
#'   `hit` (the weighted mean chance of guilt among them; `NaN` when nobody is
#'   searched).
#' @examples
#' outtst(q = c(0.1, 0.5, 0.9), t = 0.3)
#' @export
outtst <- function(q, t, w = NULL) {
  q <- .num(q, "q", lo = 0)
  if (any(q > 1)) stop("`q` must be at most 1.", call. = FALSE)
  t <- .num(t, "t", len = 1L)
  w <- if (is.null(w)) {
    rep(1, length(q))
  } else {
    .num(w, "w", len = length(q), lo = 0)
  }
  if (sum(w) <= 0) stop("`w` must have a positive total.", call. = FALSE)
  r <- .Call(C_csls_outtst, q, w, t)
  c(stop = r[1L], hit = r[2L])
}

#' Benchmark decomposition of a disparity
#'
#' Splits each group's disparity in force per resident, against a reference
#' group, into a disparity in contact per resident times a disparity in force
#' per contact. The split is a product, never a sum. See `proof("dspbm")`.
#'
#' @param pop Residents by group: positive numbers, optionally named.
#' @param contact Police contacts by group: positive numbers, same length.
#' @param force Uses of force by group: positive numbers, same length.
#' @param ref The reference group, by name or position.
#' @return A data frame with one row per group: `group`, `cd` (contact per
#'   resident), `fd` (force per contact) and `rd` (force per resident), each
#'   against the reference group, with `rd = cd * fd`.
#' @examples
#' dspbm(pop = c(a = 9000, b = 1000), contact = c(a = 900, b = 300),
#'       force = c(a = 9, b = 6))
#' @export
dspbm <- function(pop, contact, force, ref = 1) {
  nms <- .grp(pop)
  pop <- .num(pop, "pop", lo = 0, open = TRUE)
  n <- length(pop)
  contact <- .num(contact, "contact", len = n, lo = 0, open = TRUE)
  force <- .num(force, "force", len = n, lo = 0, open = TRUE)
  i <- .ref(ref, nms, n)
  r <- .Call(C_csls_dspbm, pop, contact, force, i - 1L)
  data.frame(group = nms, cd = r[["cd"]], fd = r[["fd"]], rd = r[["rd"]])
}

#' Dyadic offending rates
#'
#' Offences by one group against another, as a rate per member of the
#' offending group and as a rate per offender-victim pair exposed. The rate
#' per offender rises with the size of the victim group even when offenders
#' choose victims at random; the rate per pair does not. See
#' `proof("dyadrt")`.
#'
#' @param count A square matrix of offence counts, offending group in rows
#'   and victim group in columns: non-negative numbers.
#' @param pop Group sizes in the order of the rows: positive numbers,
#'   optionally named.
#' @return A data frame with one row per dyad: `off`, `vic`, `count`,
#'   `rate_off`, `rate_pair`, `null_off` (the rate per offender under random
#'   mixing) and `ratio` (`rate_pair` over the pooled pair rate; 1 under
#'   random mixing).
#' @examples
#' pop <- c(a = 800, b = 200)
#' dyadrt(count = outer(pop, pop) / 1000, pop = pop)
#' @export
dyadrt <- function(count, pop) {
  nms <- .grp(pop)
  pop <- .num(pop, "pop", lo = 0, open = TRUE)
  g <- length(pop)
  if (!is.matrix(count) || nrow(count) != g || ncol(count) != g) {
    stop("`count` must be a ", g, " by ", g, " matrix.", call. = FALSE)
  }
  x <- .num(as.vector(count), "count", lo = 0)
  if (sum(x) <= 0) stop("`count` holds no offences.", call. = FALSE)
  off <- rep(seq_len(g), times = g)
  vic <- rep(seq_len(g), each = g)
  r <- .Call(C_csls_dyadrt, x, off - 1L, vic - 1L, pop)
  data.frame(off = nms[off], vic = nms[vic], count = x,
             rate_off = r[["rate_off"]], rate_pair = r[["rate_pair"]],
             null_off = r[["null_off"]], ratio = r[["ratio"]])
}
