#' Clearance-rate decomposition and Simpson's paradox
#'
#' Splits the change in an aggregate clearance rate between two periods,
#' areas or services (A and B) into a rate effect, from each offence type
#' clearing more or less, and a mix effect, from the workload moving between
#' types that clear at different rates (Kitagawa's decomposition). The two
#' always add up to the change. When every type clears better in B and the
#' aggregate still falls, the fall is all mix: Simpson's paradox, and the
#' reason an aggregate clearance rate cannot be read as performance without
#' the split by type. See `proof("clrsmp")`.
#'
#' Each part's mix term is measured from the mean of the two aggregate
#' rates, so a part's term is negative when work moved toward a part that
#' clears below average; the centre changes how the mix effect is shared
#' among the parts, never its total.
#'
#' @param clra Offences cleared in A, by type: non-negative numbers,
#'   optionally named.
#' @param reca Offences recorded in A, by type: positive numbers, same length.
#' @param clrb Offences cleared in B, by type: non-negative numbers, same
#'   length.
#' @param recb Offences recorded in B, by type: positive numbers, same length.
#' @return A list. `parts`, a data frame with one row per type: `part`,
#'   `ratea`, `rateb` (clearance rates), `sharea`, `shareb` (shares of the
#'   recorded workload), `rateeff` and `mixeff`. `total`, a named numeric
#'   vector: `ratea`, `rateb` (aggregate rates), `change`, `rateeff`,
#'   `mixeff` (with `change = rateeff + mixeff`), `up` and `down` (types
#'   whose rate rose or fell) and `simpson`: 2 when every type moved one way
#'   and the aggregate the other, 1 when the aggregate change and the rate
#'   effect differ in sign, 0 otherwise.
#' @references Kitagawa, E. M. (1955). Components of a difference between
#'   two rates. *Journal of the American Statistical Association*, 50(272),
#'   1168-1194. \doi{10.1080/01621459.1955.10501299}
#' @examples
#' # both types clear better in B, yet the aggregate falls
#' clrsmp(clra = c(minor = 90, serious = 20), reca = c(100, 50),
#'        clrb = c(19, 60), recb = c(20, 140))
#' @export
clrsmp <- function(clra, reca, clrb, recb) {
  nms <- .grp(clra)
  clra <- .num(clra, "clra", lo = 0)
  n <- length(clra)
  reca <- .num(reca, "reca", len = n, lo = 0, open = TRUE)
  clrb <- .num(clrb, "clrb", len = n, lo = 0)
  recb <- .num(recb, "recb", len = n, lo = 0, open = TRUE)
  r <- .Call(C_csls_clrsmp, clra, reca, clrb, recb)
  tot <- r[["tot"]]
  cnt <- r[["cnt"]]
  list(
    parts = data.frame(part = nms, ratea = r[["ra"]], rateb = r[["rb"]],
                       sharea = r[["wa"]], shareb = r[["wb"]],
                       rateeff = r[["reff"]], mixeff = r[["meff"]]),
    total = c(ratea = tot[1L], rateb = tot[2L], change = tot[3L],
              rateeff = tot[4L], mixeff = tot[5L], up = cnt[1L],
              down = cnt[2L], simpson = cnt[3L])
  )
}

#' Rates and disparities when exposure is measured by a proxy
#'
#' Recorded counts (stops, searches, uses of force) are rates per unit of
#' exposure to police, and exposure is seldom observed. If each group's true
#' exposure lies within a factor `gamma` of an observed proxy (residents,
#' drivers on the road, calls for service), each group's true rate lies
#' within a factor `gamma` of its proxy rate, and the ratio of two groups'
#' true rates within a factor `gamma^2` of the proxy ratio. The direction of
#' a disparity is settled only when the proxy ratio lies outside that band.
#' See `proof("expbnd")`.
#'
#' @param y Recorded counts by group: positive numbers, optionally named.
#' @param m The exposure proxy by group: positive numbers, same length.
#' @param gamma The factor within which true exposure lies around the
#'   proxy: one number, at least 1.
#' @param ref The reference group, by name or position.
#' @return A data frame with one row per group: `group`, `rate` (`y / m`),
#'   `lo` and `hi` (the bounds on the true rate), `ratio` (the proxy rate
#'   against the reference group), `rlo` and `rhi` (the bounds on the true
#'   ratio; every end is reached by some admissible exposure), and `sign`: 1
#'   when the true ratio exceeds 1 for every admissible exposure, -1 when it
#'   is below 1 for every one, 0 when the data leave it open.
#' @references Knox, D., Lowe, W. and Mummolo, J. (2020). Administrative
#'   records mask racially biased policing. *American Political Science
#'   Review*, 114(3), 619-637. \doi{10.1017/S0003055420000039}
#' @examples
#' expbnd(y = c(white = 100, black = 300), m = c(1000, 1000), gamma = 1.5)
#' @export
expbnd <- function(y, m, gamma, ref = 1) {
  nms <- .grp(y)
  y <- .num(y, "y", lo = 0, open = TRUE)
  n <- length(y)
  m <- .num(m, "m", len = n, lo = 0, open = TRUE)
  gamma <- .num(gamma, "gamma", len = 1L, lo = 1)
  i <- .ref(ref, nms, n)
  r <- .Call(C_csls_expbnd, y, m, i - 1L, gamma)
  data.frame(group = nms, rate = r[["rate"]], lo = r[["lo"]], hi = r[["hi"]],
             ratio = r[["ratio"]], rlo = r[["rlo"]], rhi = r[["rhi"]],
             sign = r[["sign"]])
}

#' Risk difference when group membership is misrecorded
#'
#' When the group a person is recorded in (an officer's perception of race,
#' say) differs at random from the true one, with the same error rates
#' whatever the outcome, the measured groups are mixtures of the true ones.
#' The measured risk difference is then the true one times the factor
#' `(se + sp - 1) p (1 - p) / (d1 d0)`, where `d1` and `d0` are the shares
#' recorded in each group: with `se + sp >= 1` it never has the opposite
#' sign to the true difference and is never larger in size, so a
#' disparity in recorded data understates the true one. See `proof("misrd")`.
#'
#' @param r1 The true risk in group 1: one number from 0 to 1.
#' @param r0 The true risk in group 0: one number from 0 to 1.
#' @param prev The true share in group 1, strictly between 0 and 1.
#' @param se Sensitivity: the chance a member of group 1 is recorded in
#'   group 1, from 0 to 1.
#' @param sp Specificity: the chance a member of group 0 is recorded in
#'   group 0, from 0 to 1.
#' @return A named numeric vector: `true` (`r1 - r0`), `measured` (the
#'   difference between the recorded groups), `factor` (`measured / true`
#'   whenever `true` is not 0), `d1` and `d0`.
#' @references Bross, I. (1954). Misclassification in 2 x 2 tables.
#'   *Biometrics*, 10(4), 478-486. \doi{10.2307/3001619}
#' @examples
#' misrd(r1 = 0.3, r0 = 0.1, prev = 0.4, se = 0.8, sp = 0.9)
#' @export
misrd <- function(r1, r0, prev, se, sp) {
  r1 <- .p01(r1, "r1")
  r0 <- .p01(r0, "r0")
  prev <- .one(prev, "prev", 0, 1)
  se <- .p01(se, "se")
  sp <- .p01(sp, "sp")
  d1 <- prev * se + (1 - prev) * (1 - sp)
  if (d1 <= 0 || d1 >= 1) {
    stop("`se` and `sp` leave a recorded group empty.", call. = FALSE)
  }
  r <- .Call(C_csls_misrd, r1, r0, prev, se, sp)
  c(true = r[1L], measured = r[2L], factor = r[3L], d1 = r[4L], d0 = r[5L])
}

#' Regression to the mean at selected hot spots
#'
#' Places chosen because their count in a first period exceeds a cut show a
#' fall in the next period even with no treatment and no trend: if the two
#' periods are exchangeable, the selected places' mean change is never
#' positive, and places chosen for being low never fall on average. This
#' reports the change at the selected places next to the same statistic on
#' the data together with its period-swapped copy, which is exchangeable by
#' construction and so never positive; the excess of the first over the
#' second is the part of the fall that selection alone does not explain. A
#' treatment effect needs a control group selected the same way. See
#' `proof("rtmhs")`.
#'
#' @param x1 Counts by place in the first period: finite numbers.
#' @param x2 Counts by place in the second period: same length.
#' @param cut The selection cut: places with `x1 > cut` are selected.
#' @param w `NULL` for equal weights, or one non-negative weight per place.
#' @return A named numeric vector: `nsel` (the weight selected), `change`
#'   (the weighted mean of `x2 - x1` at the selected places), `mirror` (the
#'   mean of `x1 - x2` at places with `x2 > cut`), `sym` (the change in the
#'   data with its swapped copy; never positive), `excess` (`change - sym`),
#'   `low` (the mean change at places with `x1 < cut`) and `lowsym` (its
#'   swapped-copy version; never negative). `NaN` where no place qualifies.
#' @references Barnett, A. G., van der Pols, J. C. and Dobson, A. J. (2005).
#'   Regression to the mean: what it is and how to deal with it.
#'   *International Journal of Epidemiology*, 34(1), 215-220.
#'   \doi{10.1093/ije/dyh299}
#' @examples
#' set.seed(5)
#' mu <- rgamma(300, 2, 0.5)
#' x1 <- rpois(300, mu)
#' x2 <- rpois(300, mu)
#' rtmhs(x1, x2, cut = quantile(x1, 0.9))
#' @export
rtmhs <- function(x1, x2, cut, w = NULL) {
  x1 <- .num(x1, "x1")
  n <- length(x1)
  x2 <- .num(x2, "x2", len = n)
  cut <- .num(cut, "cut", len = 1L)
  w <- if (is.null(w)) rep(1, n) else .num(w, "w", len = n, lo = 0)
  if (sum(w[x1 > cut]) <= 0) {
    stop("no place with weight exceeds `cut`.", call. = FALSE)
  }
  r <- .Call(C_csls_rtmhs, x1, x2, w, cut)
  c(nsel = r[1L], change = r[2L], mirror = r[3L], sym = r[4L],
    excess = r[5L], low = r[6L], lowsym = r[7L])
}

# One number in the closed unit interval.
.p01 <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop("`", nm, "` must be one number from 0 to 1.", call. = FALSE)
  }
  as.double(x)
}
