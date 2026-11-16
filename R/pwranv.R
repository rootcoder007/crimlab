#' Power of the one-way analysis of variance
#'
#' The power of the F test of equal means in `k` groups of `n` with
#' Cohen's effect size `f` (the standard deviation of the group means
#' about their mean, over the common standard deviation), from the
#' noncentral F distribution; or, given the power wanted, the smallest
#' whole `n` per group that reaches it.
#'
#' The statistic has `df1 = k - 1` and `df2 = k (n - 1)` degrees of
#' freedom and noncentrality `ncp = k n f^2` (Cohen 1988, ch. 8); the
#' power is the chance that it passes `fcrit`, the upper `alpha` quantile
#' of the central F. The noncentral upper tail is the Poisson mixture of
#' beta upper tails (Johnson, Kotz and Balakrishnan 1995, ch. 30). With
#' `power`, `n` is searched over the whole numbers from 2 to 10000000;
#' where no such `n` reaches it (`f = 0`, or too small an effect) the row
#' is `NA` from `df2` on. Cohen's conventions are 0.1, 0.25 and 0.4 for a
#' small, medium and large effect.
#'
#' @param k The number of groups: a whole number from 2 to 1000000.
#' @param n `NULL`, or the sizes of each group to try: numbers above 1, at
#'   most 10000000.
#' @param f The effect sizes to try: numbers from 0.
#' @inheritParams pwrsim
#' @param power `NULL`, or the power wanted: one number above `alpha` and
#'   below 1. Give one of `n` and `power`.
#' @return A data frame, one row per `n` and `f` (`n` varying fastest), or
#'   per `f` with `power`: `k`, `n`, `f`, `df1`, `df2`, `ncp`, `fcrit` and
#'   `power`.
#' @references
#'   Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
#'   Sciences*, 2nd ed. Lawrence Erlbaum. \doi{10.4324/9780203771587}
#'
#'   Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). *Continuous
#'   Univariate Distributions*, vol. 2, 2nd ed. Wiley.
#' @examples
#' pwranv(3, n = c(20, 40, 60), f = 0.25)
#' pwranv(4, f = c(0.1, 0.25, 0.4), power = 0.8)
#' @export
pwranv <- function(k, n = NULL, f, alpha = 0.05, power = NULL) {
  fn <- "pwranv"
  k <- .cknum(k, "k", fn, min = 2, max = 1e6, whole = TRUE)
  fv <- .cknum(f, "f", fn, scalar = FALSE, min = 0)
  alpha <- .cknum(alpha, "alpha", fn)
  if (alpha <= 0 || alpha >= 1) {
    .ckstop(fn, "`alpha` must be one number strictly between 0 and 1")
  }
  if (is.null(n) == is.null(power)) {
    .ckstop(fn, "give one of `n` and `power`")
  }
  if (!is.null(n)) {
    nv <- .cknum(n, "n", fn, scalar = FALSE, min = 1, max = 1e7)
    if (any(nv <= 1)) .ckstop(fn, "`n` must be above 1")
    g <- expand.grid(n = nv, f = fv)
    nv <- g$n
    fv <- g$f
    tg <- NA_real_
  } else {
    tg <- .cknum(power, "power", fn)
    if (tg <= alpha || tg >= 1) {
      .ckstop(fn, "`power` must be above `alpha` and below 1")
    }
    nv <- NULL
  }
  o <- .Call(C_clab_r_pwranv, k, nv, fv, alpha, tg)
  m <- length(fv)
  o[is.nan(o)] <- NA_real_
  data.frame(k = rep(k, m), n = o[seq_len(m)], f = fv,
             df1 = o[m + seq_len(m)], df2 = o[2 * m + seq_len(m)],
             ncp = o[3 * m + seq_len(m)], fcrit = o[4 * m + seq_len(m)],
             power = o[5 * m + seq_len(m)])
}
