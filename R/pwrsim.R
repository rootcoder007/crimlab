#' Power by simulation
#'
#' The power of a two-arm trial found by simulating it: for each control
#' arm size in `n` and each effect in `eff`, `nsim` trials are drawn, each
#' is tested, and the power is the share rejected at `alpha`, with its
#' Monte Carlo standard error and Wilson's interval (Morris, White and
#' Crowther 2019). It answers what a formula cannot: an outcome that is
#' binary or a count, unequal arms, a permutation test.
#'
#' The control arm has `n` units and the treated arm `round(ratio * n)`.
#' The outcome is normal with standard deviation `sd` and mean 0 in the
#' control arm and `eff` in the treated arm (`fam = "normal"`), 0 or 1
#' with chance `p0` and `p0 + eff` (`"binary"`), or Poisson with mean `mu0`
#' and `mu0 + eff` (`"poisson"`), drawn by inversion from the package's
#' generator. The test is Welch's t test (`test = "welch"`; with no
#' spread in either arm the p-value is 0 when the means differ in the
#' direction tested, else 1) or the permutation test of the difference in
#' means with `nperm` relabellings (`"perm"`; `p = (1 + k) / (1 + nperm)`).
#' Each design has a stream of its own made from `seed` and its place in
#' the grid, so a design's power does not depend on the others.
#'
#' @param n The control arm sizes to try: whole numbers from 2 to 1000000.
#' @param eff The effects to try: differences in means (in the outcome's
#'   units; risk differences for `"binary"`).
#' @param fam The outcome: `"normal"`, `"binary"` or `"poisson"`.
#' @param sd The standard deviation of a normal outcome: positive.
#' @param p0 The control arm's chance of a binary outcome, from 0 to 1;
#'   `p0 + eff` must be too.
#' @param mu0 The control arm's mean of a Poisson outcome: positive;
#'   `mu0 + eff` must be at least 0.
#' @param ratio The treated arm's size over the control arm's: positive.
#' @param test `"welch"` or `"perm"`.
#' @param alt The alternative: `"two.sided"`, `"greater"` (the treated
#'   higher) or `"less"`.
#' @param alpha The level of the test: one number strictly between 0 and
#'   1.
#' @param nsim The trials drawn per design: a whole number from 1 to
#'   1000000.
#' @param nperm The relabellings of the permutation test: a whole number
#'   from 1 to 100000 with `1 / (nperm + 1)` at most `alpha`.
#' @param level The level of Wilson's interval: one number strictly
#'   between 0 and 1.
#' @inheritParams rndalc
#' @return An object of class `pwrsim`, a list: `tab`, a data frame with
#'   one row per design (`n` varying fastest): `n`, `n1` (the treated
#'   arm), `eff`, `power`, `se`, `lo`, `hi` and `rej` (the trials
#'   rejected); and `fam`, `test`, `alt`, `alpha`, `nsim` and `seed`.
#' @references
#'   Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using
#'   simulation studies to evaluate statistical methods. *Statistics in
#'   Medicine*, 38(11), 2074-2102. \doi{10.1002/sim.8086}
#'
#'   Welch, B. L. (1947). The generalization of 'Student's' problem when
#'   several different population variances are involved. *Biometrika*,
#'   34(1-2), 28-35. \doi{10.1093/biomet/34.1-2.28}
#'
#'   Wilson, E. B. (1927). Probable inference, the law of succession, and
#'   statistical inference. *Journal of the American Statistical
#'   Association*, 22(158), 209-212. \doi{10.1080/01621459.1927.10502953}
#' @examples
#' pwrsim(c(20, 50), eff = c(0, 0.5), nsim = 400, seed = 1)
#' pwrsim(100, eff = 0.15, fam = "binary", p0 = 0.3, nsim = 400,
#'        seed = 2)$tab
#' @export
pwrsim <- function(
  n,
  eff,
  fam = c("normal", "binary", "poisson"),
  sd = 1,
  p0 = 0.5,
  mu0 = 1,
  ratio = 1,
  test = c("welch", "perm"),
  alt = c("two.sided", "greater", "less"),
  alpha = 0.05,
  nsim = 1000,
  nperm = 199,
  level = 0.95,
  seed = NULL
) {
  f <- "pwrsim"
  n <- .cknum(n, "n", f, scalar = FALSE, min = 2, max = 1e6, whole = TRUE)
  eff <- .cknum(eff, "eff", f, scalar = FALSE)
  fam <- .ckmarg("fam", f)
  test <- .ckmarg("test", f)
  alt <- .ckmarg("alt", f)
  par <- switch(
    fam,
    normal = .cknum(sd, "sd", f, positive = TRUE),
    binary = .cknum(p0, "p0", f, min = 0, max = 1),
    poisson = .cknum(mu0, "mu0", f, positive = TRUE)
  )
  if (fam == "binary" && any(par + eff < 0 | par + eff > 1)) {
    .ckstop(f, "`p0 + eff` must be from 0 to 1")
  }
  if (fam == "poisson" && any(par + eff < 0)) {
    .ckstop(f, "`mu0 + eff` must be at least 0")
  }
  ratio <- .cknum(ratio, "ratio", f, positive = TRUE)
  alpha <- .xpopen(alpha, "alpha", f)
  nsim <- .cknum(nsim, "nsim", f, min = 1, max = 1e6, whole = TRUE)
  if (test == "perm") {
    nperm <- .cknum(nperm, "nperm", f, min = 1, max = 1e5, whole = TRUE)
    if (1 / (nperm + 1) > alpha) {
      .ckstop(f, paste("`nperm` is too few for `alpha`: the least p-value,",
                       "1 / (nperm + 1), is above it"))
    }
  }
  level <- .xpopen(level, "level", f)
  seed <- .xpseed(seed, f)
  g <- expand.grid(n = n, eff = eff)
  n1 <- round(ratio * g$n)
  if (any(n1 < 2 | n1 > 1e6)) {
    .ckstop(f, "`ratio * n` must round to treated arms of 2 to 1000000")
  }
  fk <- match(fam, c("normal", "binary", "poisson")) - 1L
  o <- .Call(C_clab_r_pwrsim, fk, par, as.integer(g$n), as.integer(n1), g$eff,
             match(test, c("welch", "perm")) - 1L,
             match(alt, c("two.sided", "greater", "less")) - 1L, alpha,
             as.integer(nsim), as.integer(if (test == "perm") nperm else 1),
             level, seed)
  m <- nrow(g)
  tab <- data.frame(n = g$n, n1 = n1, eff = g$eff, power = o[seq_len(m)],
                    se = o[m + seq_len(m)], lo = o[2 * m + seq_len(m)],
                    hi = o[3 * m + seq_len(m)], rej = o[4 * m + seq_len(m)])
  structure(list(tab = tab, fam = fam, test = test, alt = alt, alpha = alpha,
                 nsim = nsim, seed = seed), class = "pwrsim")
}

#' @rdname pwrsim
#' @param x An object of class `pwrsim`.
#' @inheritParams cemtch
#' @export
print.pwrsim <- function(x, ...) {
  cat("pwrsim: ", x[["fam"]], " outcome, ", x[["test"]], " test (",
      x[["alt"]], ") at alpha ", format(x[["alpha"]]), ", ",
      format(x[["nsim"]], scientific = FALSE), " trials per design\n",
      sep = "")
  print(x[["tab"]], row.names = FALSE)
  invisible(x)
}

#' One number strictly between 0 and 1.
#' @noRd
.xpopen <- function(x, what, f) {
  x <- .cknum(x, what, f)
  if (x <= 0 || x >= 1) {
    .ckstop(f, "`%s` must be one number strictly between 0 and 1", what)
  }
  x
}
