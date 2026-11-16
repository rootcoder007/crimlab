#' Permutation test within blocks
#'
#' The randomisation test of a treatment effect in a blocked (or
#' stratified, or matched) experiment: the treatment labels are permuted
#' within each block, so that every arrangement keeps each block's number
#' of treated units, as the design randomised them, and the observed
#' statistic is compared with its values over the arrangements (Rosenbaum
#' 2002, ch. 2; Good 2005, sec. 3.6).
#'
#' `stat = "diff"` is the difference between the mean outcome of the
#' treated and of the controls over all units; `stat = "strat"` is the
#' mean over blocks of the blocks' own differences, weighted by the
#' blocks' sizes, the stratified estimate of the average effect (every
#' block then needs treated and control units). With `exact = TRUE` every
#' arrangement is enumerated (at most 1000000) and `p` is the share at
#' least as extreme as the observed value; otherwise `nperm` arrangements
#' are drawn from the package's generator and `p = (1 + k) / (1 + nperm)`,
#' `k` of them at least as extreme, which is never 0 (Phipson and Smyth
#' 2010). Values within `1e-9 (|T| + max|y|)` of the observed `T` count
#' as at least as extreme, so that arrangements equal to it in exact
#' arithmetic are never lost to rounding.
#'
#' @inheritParams cemtch
#' @inheritParams ifesc
#' @inheritParams pwrsim
#' @param d The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`),
#'   both present.
#' @param block The name of the block column (labels).
#' @param stat `"diff"` (the difference in means) or `"strat"` (the
#'   block-weighted mean of the blocks' differences).
#' @param nperm The number of arrangements drawn: a whole number from 1 to
#'   10000000 (not used with `exact`).
#' @param exact `TRUE` to enumerate every arrangement.
#' @param seed `NULL` (drawn from R's generator, so [set.seed()] fixes it)
#'   or a whole number from 0 to 2147483647; not used with `exact`.
#' @return An object of class `blkprm`, a list: `est` (the observed
#'   statistic), `p`, `nperm` (the arrangements used), `nge` (those at
#'   least as extreme), `nblock`, `stat`, `alt`, `exact`, `dist` (the
#'   permuted values) and `seed` (`NA` with `exact`).
#' @references
#'   Rosenbaum, P. R. (2002). *Observational Studies*, 2nd ed. Springer.
#'   \doi{10.1007/978-1-4757-3692-2}
#'
#'   Good, P. I. (2005). *Permutation, Parametric, and Bootstrap Tests of
#'   Hypotheses*, 3rd ed. Springer. \doi{10.1007/b138696}
#'
#'   Phipson, B. and Smyth, G. K. (2010). Permutation p-values should never
#'   be zero: calculating exact p-values when permutations are randomly
#'   drawn. *Statistical Applications in Genetics and Molecular Biology*,
#'   9(1), 39. \doi{10.2202/1544-6115.1585}
#' @examples
#' set.seed(4)
#' d <- data.frame(site = rep(1:6, each = 4), t = rep(c(1, 0, 1, 0), 6))
#' d$y <- d$site + 0.8 * d$t + rnorm(24, sd = 0.5)
#' blkprm(d, "y", "t", "site", nperm = 999, seed = 1)
#' blkprm(d[1:12, ], "y", "t", "site", stat = "strat", exact = TRUE)$p
#' @export
blkprm <- function(
  data,
  y,
  d,
  block,
  stat = c("diff", "strat"),
  alt = c("two.sided", "greater", "less"),
  nperm = 2000,
  exact = FALSE,
  seed = NULL
) {
  f <- "blkprm"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f, n = 1L)
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, block, "block", f, n = 1L)
  stat <- .ckmarg("stat", f)
  alt <- .ckmarg("alt", f)
  yv <- .axmat(dd, y, f)[, 1L]
  tr <- .axmat(dd, d, f)[, 1L]
  if (!all(tr %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", d)
  if (all(tr == 1) || all(tr == 0)) {
    .ckstop(f, "column `%s` must hold treated and control units", d)
  }
  b <- .axcode(dd[[block]])
  m <- tabulate(b$c + 1L, b$n)
  t1 <- tabulate(b$c[tr == 1] + 1L, b$n)
  if (stat == "strat" && any(t1 == 0L | t1 == m)) {
    .ckstop(f, paste("with `stat = \"strat\"` every block needs treated and",
                     "control units"))
  }
  exact <- .ckflag(exact, "exact", f)
  if (exact) {
    tot <- prod(choose(m, t1))
    if (tot > 1e6) {
      .ckstop(f, paste("`exact` enumerates at most 1000000 arrangements;",
                       "these blocks have %s"), format(tot))
    }
    np <- tot
    seed <- NA_real_
  } else {
    np <- .cknum(nperm, "nperm", f, min = 1, max = 1e7, whole = TRUE)
    seed <- .xpseed(seed, f)
  }
  k <- .Call(C_clab_r_blkprm, yv, as.integer(tr), as.integer(b$c), b$n,
             match(stat, c("diff", "strat")) - 1L,
             match(alt, c("two.sided", "greater", "less")) - 1L,
             as.integer(np), exact, seed)
  st <- k$st
  structure(
    list(est = st[1L], p = st[2L], nperm = st[3L], nge = st[4L],
         nblock = b$n, stat = stat, alt = alt, exact = exact, dist = k$dist,
         seed = seed),
    class = "blkprm"
  )
}

#' @rdname blkprm
#' @param x An object of class `blkprm`.
#' @inheritParams cemtch
#' @export
print.blkprm <- function(x, ...) {
  cat("blkprm: ", if (x[["stat"]] == "diff") "difference in means" else
        "block-weighted difference", " ", format(x[["est"]], digits = 4),
      ", p = ", format(x[["p"]], digits = 3), " (", x[["alt"]], ", ",
      if (x[["exact"]]) "all " else "",
      format(x[["nperm"]], scientific = FALSE),
      " arrangements within ", x[["nblock"]], " blocks)\n", sep = "")
  invisible(x)
}
