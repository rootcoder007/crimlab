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
#' @param day Stops by group in daylight: non-negative numbers, as a numeric
#'   vector, a 1-d table, a list of single numbers or a one-row data frame.
#'   With `night` left out, one table of groups by stops: a matrix, 2-d table
#'   or data frame whose two columns are named `day` and `night`, in either
#'   order (columns are never taken by position).
#' @param night Stops by group in darkness: positive numbers, in the forms of
#'   `day`; named inputs are matched to `day` by name, and an error lists any
#'   group named in one and not the other. Unnamed inputs are matched by
#'   position: only names protect against a mismatch.
#' @param ref The reference group, by name (matched exactly, case included;
#'   a factor by its label) or by position in `day` as given.
#' @param shift `NULL`, or one positive number per group: daylight exposure
#'   over darkness exposure, in the forms of `day`, matched by name.
#' @return A data frame with one row per group, in the order of `day`:
#'   `group`, `ref` (`TRUE` for the reference group), `shday`,
#'   `shnight`, `ratio`, `sign` (the sign of `shday - shnight`)
#'   and, when `shift` is given, `adj`.
#' @references Grogger, J. and Ridgeway, G. (2006). Testing for racial
#'   profiling in traffic stops from behind a veil of darkness. *Journal of
#'   the American Statistical Association*, 101(475), 878-887.
#'   \doi{10.1198/016214506000000168}
#' @examples
#' vodtst(day = c(white = 600, black = 400),
#'        night = c(white = 700, black = 300))
#' # one table of groups by day and night
#' stops <- data.frame(day = c(600, 400), night = c(700, 300),
#'                     row.names = c("white", "black"))
#' vodtst(stops)
#' @export
vodtst <- function(day, night = NULL, ref = 1, shift = NULL) {
  .Call(C_clab_r_vodtst, day, night, ref, shift)
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
  .Call(C_clab_r_orrisk, or, rate, share)
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
#'   of `q` with a positive total; when both are named with unique names,
#'   `w` is matched to `q` by name.
#' @return A named numeric vector: `stop` (the weighted share searched) and
#'   `hit` (the weighted mean chance of guilt among them; `NaN` when nobody is
#'   searched).
#' @examples
#' outtst(q = c(0.1, 0.5, 0.9), t = 0.3)
#' @export
outtst <- function(q, t, w = NULL) {
  .Call(C_clab_r_outtst, q, t, w)
}

#' Benchmark decomposition of a disparity
#'
#' Splits each group's disparity in force per resident, against a reference
#' group, into a disparity in contact per resident times a disparity in force
#' per contact. The split is a product, never a sum. See `proof("dspbm")`.
#'
#' @param pop Residents by group: positive numbers, as a numeric vector, a
#'   1-d table, a list of single numbers or a one-row data frame.
#' @param contact Police contacts by group: positive numbers, in the forms of
#'   `pop`; named inputs are matched to `pop` by name, unnamed ones by
#'   position (only names protect against a mismatch).
#' @param force Uses of force by group: positive numbers, in the forms of
#'   `pop`, matched by name.
#' @param ref The reference group, by name (matched exactly, case included;
#'   a factor by its label) or by position in `pop` as given.
#' @return A data frame with one row per group, in the order of `pop`:
#'   `group`, `ref` (`TRUE` for the reference group), `cd` (contact per
#'   resident), `fd` (force per contact) and `rd` (force per resident), each
#'   against the reference group, with `rd = cd * fd`.
#' @examples
#' dspbm(pop = c(a = 9000, b = 1000), contact = c(a = 900, b = 300),
#'       force = c(a = 9, b = 6))
#' @export
dspbm <- function(pop, contact, force, ref = 1) {
  .Call(C_clab_r_dspbm, pop, contact, force, ref)
}

#' Dyadic offending rates
#'
#' Offences by one group against another, as a rate per member of the
#' offending group and as a rate per offender-victim pair exposed. The rate
#' per offender rises with the size of the victim group even when offenders
#' choose victims at random; the rate per pair does not. See
#' `proof("dyadrt")`.
#'
#' @param count A square matrix (or 2-d table) of offence counts, offending
#'   group in rows and victim group in columns: non-negative numbers. When
#'   it names both its rows and its columns, the columns are matched to the
#'   rows by name, and both to `pop` when `pop` is named.
#' @param pop Group sizes: positive numbers, as a numeric vector, a 1-d
#'   table, a list of single numbers or a one-row data frame. Unnamed, they
#'   are taken in the order of the rows: only names protect against a
#'   mismatch.
#' @return A data frame with one row per dyad: `off`, `vic`, `count`,
#'   `rtoff`, `rtpair`, `nulloff` (the rate per offender under random
#'   mixing) and `ratio` (`rtpair` over the pooled pair rate; 1 under
#'   random mixing).
#' @examples
#' pop <- c(a = 800, b = 200)
#' dyadrt(count = outer(pop, pop) / 1000, pop = pop)
#' @export
dyadrt <- function(count, pop) {
  .Call(C_clab_r_dyadrt, count, pop)
}
