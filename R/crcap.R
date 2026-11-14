#' Capture-recapture estimates of a hidden population
#'
#' How many offenders, victims or incidents no list records, from the
#' overlap of two or three lists (police records, hospital records, a
#' survey). With two lists the Lincoln-Petersen estimate `n1 n2 / m`
#' assumes the lists are independent; when being on one list changes the
#' odds of being on the other by at most a factor `kappa`, the truth lies
#' within `[n1 n2 / (kappa m), kappa n1 n2 / m]` (`lower`, `upper`), and no
#' data can narrow that box. Chapman's estimate `(n1 + 1)(n2 + 1) / (m + 1)
#' - 1` is nearly unbiased, with Seber's standard error.
#'
#' With three lists the seven observed cells identify the never-listed cell
#' only under an assumption on the three-way interaction of a log-linear
#' model: with none, it is `m100 m010 m001 m111 / (m110 m101 m011)`, and a
#' three-way log interaction `t` multiplies it by `exp(-t)`; `int` gives
#' that sensitivity. The pairwise estimates are given beside it, each at
#' least the number seen on either of its two lists (`floor`). See
#' `proof("crcap")`.
#'
#' @param x `c(n1, n2, m)`, the sizes of two lists and their overlap
#'   (named so or in that order), or the seven cells of three lists, named
#'   by membership: `"100"`, `"010"`, `"001"`, `"110"`, `"101"`, `"011"`,
#'   `"111"`. Counts from 0.
#' @param kappa The largest factor by which being on one list changes the
#'   odds of being on the other: one number from 1 (two lists).
#' @param int `NULL`, or three-way log interactions at which to give the
#'   missing cell (three lists).
#' @return An object of class `crcap`, a list: `type` (`"two"` or
#'   `"three"`); `est`, a named vector (two lists: `petersen`, `lower`,
#'   `upper`, `chapman`, `se`, `floor`; three: `observed`, `missing`,
#'   `total`); `pairs` (three lists: `lists`, `n1`, `n2`, `m`, `floor`,
#'   `petersen`, `chapman`) and `sens` (with `int`: `int`, `missing`,
#'   `total`), otherwise `NULL`.
#' @references
#'   Chapman, D. G. (1951). Some properties of the hypergeometric
#'   distribution with applications to zoological sample censuses.
#'   *University of California Publications in Statistics*, 1, 131-160.
#'
#'   Seber, G. A. F. (1970). The effects of trap response on tag recapture
#'   estimates. *Biometrics*, 26(1), 13-22.
#'   \doi{10.2307/2529040}
#'
#'   Fienberg, S. E. (1972). The multiple recapture census for closed
#'   populations and incomplete 2^k contingency tables. *Biometrika*,
#'   59(3), 591-603. \doi{10.1093/biomet/59.3.591}
#' @examples
#' crcap(c(n1 = 500, n2 = 400, m = 80), kappa = 1.5)
#' crcap(c("100" = 120, "010" = 90, "001" = 60, "110" = 30, "101" = 20,
#'         "011" = 15, "111" = 8), int = c(-0.5, 0, 0.5))
#' @export
crcap <- function(x, kappa = 1, int = NULL) {
  structure(.Call(C_clab_r_crcap, x, kappa, int), class = "crcap")
}

#' The dark figure of crime: bounds instead of a point
#'
#' What recorded crime and a victimisation survey together say about the
#' crime that is never recorded, as intervals that hold under stated
#' limits on the errors (Manski 2003).
#'
#' * With a survey rate `v` and a recorded rate `r`, and survey false
#'   positives at most `amax` and false negatives at most `bmax`, the true
#'   rate lies in `[max(r, (v - amax) / (1 - amax)), min(1, v / (1 -
#'   bmax))]` and the dark figure (true minus recorded) in that interval
#'   shifted by `r`. Given `thr`, `brk` is the largest `bmax` under which
#'   the conclusion "the true rate is below `thr`" survives: `1 - v / thr`.
#' * When the survey covers only some offences (type A, `va` of them, `ra`
#'   recorded) and the recorded count of the rest is `rb`, recorded at no
#'   less than `lam` times A's rate, the total and the dark figure have the
#'   bounds in `cmp`.
#' * With `k`, the number of offences in each recorded incident, the
#'   hierarchy rule's arithmetic: one incident in the counts stands for
#'   between 1 and `max(k)` offences.
#'
#' See `proof("dkfig")`.
#'
#' @param v `NULL`, or survey (true) victimisation rates: numbers from 0 to
#'   1.
#' @param r The recorded rates: one number or one per `v`, from 0 to 1.
#' @param amax,bmax The largest false-positive and false-negative rates of
#'   the survey: numbers from 0, with a sum below 1.
#' @param thr `NULL`, or a threshold above every `v`.
#' @param va,ra,rb `NULL`, or the survey's offences of the covered type,
#'   the recorded ones of that type (both above 0, `ra` at most `va`) and
#'   the recorded offences of the uncovered type (from 0).
#' @param lam The least ratio of the uncovered type's recording rate to the
#'   covered type's: above 0, at most 1.
#' @param k `NULL`, or the number of offences in each incident: whole
#'   numbers from 1.
#' @return An object of class `dkfig`, a list: `bnd` (with `v`: `v`, `r`,
#'   `lo`, `hi`, `dlo`, `dhi` the dark figure's bounds, `ratio` the upper
#'   bound over the recorded rate, and with `thr` `brk`), `cmp` (with `va`:
#'   `rate` A's recording rate, `rateup`, the total's bounds `totlo`,
#'   `tothi`, the dark figure's `dklo`, `dkhi` and `width`) and `hier`
#'   (with `k`: `incid`, `offen`, `ratio`, `extra`, `lower`, `upper`); the
#'   parts not asked for are `NULL`.
#' @references Manski, C. F. (2003). *Partial Identification of Probability
#'   Distributions*. Springer. \doi{10.1007/b97478}
#' @examples
#' dkfig(v = c(0.12, 0.08), r = 0.03, amax = 0.01, bmax = 0.2, thr = 0.15)
#' dkfig(va = 900, ra = 300, rb = 450, lam = 0.5)
#' dkfig(k = c(1, 1, 2, 1, 3))
#' @export
dkfig <- function(
  v = NULL,
  r = NULL,
  amax = 0,
  bmax = 0,
  thr = NULL,
  va = NULL,
  ra = NULL,
  rb = NULL,
  lam = 1,
  k = NULL
) {
  o <- .Call(C_clab_r_dkfig, v, r, amax, bmax, thr, va, ra, rb, lam, k)
  structure(o, class = "dkfig")
}

#' @export
print.crcap <- function(x, ...) {
  cat("crcap:", if (x[["type"]] == "two") "two lists" else "three lists", "\n")
  print(x[["est"]])
  if (!is.null(x[["pairs"]])) {
    print(x[["pairs"]], row.names = FALSE)
  }
  if (!is.null(x[["sens"]])) {
    print(x[["sens"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.dkfig <- function(x, ...) {
  cat("dkfig: bounds on the dark figure\n")
  for (p in c("bnd", "cmp", "hier")) {
    if (!is.null(x[[p]])) print(x[[p]])
  }
  invisible(x)
}
