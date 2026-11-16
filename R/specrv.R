#' Specification curve
#'
#' The estimate of a treatment effect under every reasonable specification
#' at once, and whether the specifications as a whole show an effect: the
#' specification curve of Simonsohn, Simmons and Nelson (2020). A finding
#' that holds under one specification only is fragile; one that survives
#' the choices an analyst could defensibly have made is not.
#'
#' The specifications are every outcome in `y` crossed with every subset
#' of the optional covariates `opt` (`2^length(opt)` subsets, from none to
#' all), the covariates `x` always in. Each is the least squares fit of
#' the outcome on an intercept, the covariates and `d`, giving the
#' coefficient of `d` with a classical or HC1 standard error and
#' Student's t. The curve is summarised three ways (their sec. 3): the
#' median estimate; the number of specifications significant at `1 -
#' level` in the dominant direction (the sign of the median); and the sum
#' of the specifications' t statistics over the square root of their
#' number. Each summary is compared with its values on `nsim` data sets
#' where the null holds: with `null = "perm"` the treatment is permuted
#' across rows (for a randomised treatment); with `null = "boot"` each
#' specification's estimate is removed from its outcome (`y - b d`) and
#' rows are drawn with replacement. The p-value of a summary is `(1 + k) /
#' (1 + nok)`, `k` of the `nok` usable data sets at least as extreme in
#' the dominant direction; a bootstrap draw whose fit is singular is left
#' out.
#'
#' @inheritParams rndalc
#' @inheritParams plntfx
#' @param y The names of the outcome columns (numbers): each is a
#'   specification of the outcome.
#' @param x `NULL` or the names of the covariates in every specification
#'   (numbers); for `print()`, an object of class `specrv`.
#' @param opt `NULL` or the names of the optional covariates (numbers), at
#'   most ten.
#' @param level The level of the intervals and of significance: one
#'   number strictly between 0 and 1.
#' @param null `"perm"` or `"boot"`.
#' @param nsim The data sets under the null: a whole number from 1 to
#'   100000.
#' @return An object of class `specrv`, a list: `spec`, one row per
#'   specification (`outcome`, `ctrl`, the optional covariates in it or
#'   `"none"`, `est`, `se`, `stat`, `p`, `lo`, `hi`, `sig`); `test`, one
#'   row per summary (`stat`: `median`, `nsig`, `z`; `obs`; `p`); `sign`,
#'   the dominant direction; `dist`, the summaries under the null (a
#'   matrix, `NA` for a draw left out); `nok`; and `null`, `nsim`, `seed`.
#' @references
#'   Simonsohn, U., Simmons, J. P. and Nelson, L. D. (2020). Specification
#'   curve analysis. *Nature Human Behaviour*, 4(11), 1208-1214.
#'   \doi{10.1038/s41562-020-0912-z}
#'
#'   Steegen, S., Tuerlinckx, F., Gelman, A. and Vanpaemel, W. (2016).
#'   Increasing transparency through a multiverse analysis. *Perspectives
#'   on Psychological Science*, 11(5), 702-712.
#'   \doi{10.1177/1745691616658637}
#' @examples
#' set.seed(7)
#' n <- 200
#' d <- data.frame(t = rbinom(n, 1, 0.5), age = runif(n, 18, 60),
#'                 prior = rpois(n, 2), male = rbinom(n, 1, 0.5))
#' d$arrest <- 0.3 * d$t + 0.1 * d$prior + rnorm(n)
#' d$charge <- 0.2 * d$t + 0.1 * d$prior + rnorm(n)
#' s <- specrv(d, c("arrest", "charge"), "t", opt = c("age", "prior",
#'             "male"), nsim = 200, seed = 1)
#' s
#' head(s$spec)
#' @export
specrv <- function(
  data,
  y,
  d,
  x = NULL,
  opt = NULL,
  se = c("hc1", "iid"),
  level = 0.95,
  null = c("perm", "boot"),
  nsim = 500,
  seed = NULL
) {
  f <- "specrv"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f)
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, x, "x", f, null = TRUE)
  .axcols(dd, opt, "opt", f, null = TRUE)
  a <- c(y, d, x, opt)
  if (anyDuplicated(a)) {
    .ckstop(f, "`y`, `d`, `x` and `opt` must name different columns")
  }
  if (length(opt) > 10L) .ckstop(f, "`opt` may name at most ten columns")
  se <- .ckmarg("se", f)
  null <- .ckmarg("null", f)
  ym <- .axmat(dd, y, f)
  dv <- .axmat(dd, d, f)[, 1L]
  xm <- .axmat(dd, c(x, opt), f)
  n <- nrow(dd)
  if (n <= ncol(xm) + 2L) {
    .ckstop(f, "`data` needs more rows than covariates plus two")
  }
  if (all(dv == dv[1L])) .ckstop(f, "column `%s` must vary", d)
  ns <- 2^length(opt)
  if (ns * n * (ncol(xm) + 1) + ns * length(y) * n > 2^26) {
    .ckstop(f, "too many specifications for these rows")
  }
  level <- .xpopen(level, "level", f)
  nsim <- .cknum(nsim, "nsim", f, min = 1, max = 1e5, whole = TRUE)
  seed <- .xpseed(seed, f)
  k <- .Call(C_clab_r_specrv, xm, length(x), length(opt), ym, dv,
             se == "hc1", level, null == "boot", as.integer(nsim), seed)
  K <- length(y) * ns
  s <- k$spec
  ctrl <- vapply(seq_len(ns) - 1L, function(b) {
    u <- opt[bitwAnd(b, as.integer(2^(seq_along(opt) - 1L))) > 0L]
    if (length(u)) paste(u, collapse = " + ") else "none"
  }, "")
  p <- s[3L * K + seq_len(K)]
  spec <- data.frame(outcome = rep(y, each = ns), ctrl = rep(ctrl, length(y)),
                     est = s[seq_len(K)], se = s[K + seq_len(K)],
                     stat = s[2L * K + seq_len(K)], p = p,
                     lo = s[4L * K + seq_len(K)], hi = s[5L * K + seq_len(K)],
                     sig = p < 1 - level, stringsAsFactors = FALSE)
  st <- k$st
  dist <- matrix(k$dist, nsim, 3L,
                 dimnames = list(NULL, c("median", "nsig", "z")))
  dist[is.nan(dist)] <- NA_real_
  structure(
    list(spec = spec,
         test = data.frame(stat = c("median", "nsig", "z"), obs = st[1:3],
                           p = st[4:6], stringsAsFactors = FALSE),
         sign = st[7L], dist = dist, nok = st[8L], null = null, nsim = nsim,
         seed = seed),
    class = "specrv"
  )
}

#' @rdname specrv
#' @inheritParams cemtch
#' @export
print.specrv <- function(x, ...) {
  cat("specrv: ", nrow(x[["spec"]]), " specifications, ",
      sum(x[["spec"]]$sig), " significant; ", x[["null"]], " null, ",
      format(x[["nok"]], scientific = FALSE), " data sets\n", sep = "")
  print(x[["test"]], row.names = FALSE, digits = 4)
  invisible(x)
}
