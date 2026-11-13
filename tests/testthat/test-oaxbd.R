# oaxbd(): against lm() in each group and the formulas of Oaxaca, Blinder
# and Jann written here, the exact identities, the earlier package's
# recorded values, edge cases and refusals.

ref <- readRDS(test_path("ref-fair.rds"))
od <- ref$oax$d

# the decomposition and its delta-method variances, from lm() and cov()
byhand <- function(y, x, g, b) {
  a <- setdiff(unique(g), b)
  fit <- function(k) {
    lm(y[g == k] ~ ., data = as.data.frame(x)[g == k, , drop = FALSE])
  }
  fa <- fit(a)
  fb <- fit(b)
  xa <- c(1, colMeans(as.matrix(x)[g == a, , drop = FALSE]))
  xb <- c(1, colMeans(as.matrix(x)[g == b, , drop = FALSE]))
  vx <- function(k) {
    m <- matrix(0, length(xa), length(xa))
    m[-1, -1] <- cov(as.matrix(x)[g == k, , drop = FALSE]) / sum(g == k)
    m
  }
  sa <- vx(a)
  sb <- vx(b)
  ba <- coef(fa)
  bb <- coef(fb)
  va <- vcov(fa)
  vb <- vcov(fb)
  dd <- xa - xb
  dl <- ba - bb
  q <- function(u, m, w) drop(t(u) %*% m %*% w)
  list(gap = mean(y[g == a]) - mean(y[g == b]),
       exp = sum(dd * bb), unexp = sum(xa * dl), expa = sum(dd * ba),
       unexpa = sum(xb * dl), inter = sum(dd * dl),
       se = sqrt(c(exp = q(dd, vb, dd) + q(bb, sa + sb, bb),
                   unexp = q(xa, va + vb, xa) + q(dl, sa, dl),
                   expa = q(dd, va, dd) + q(ba, sa + sb, ba),
                   unexpa = q(xb, va + vb, xb) + q(dl, sb, dl),
                   inter = q(dd, va + vb, dd) + q(dl, sa + sb, dl))),
       segap = sqrt(var(y[g == a]) / sum(g == a) +
                      var(y[g == b]) / sum(g == b)),
       ba = ba, bb = bb, va = va, vb = vb, xa = xa, xb = xb, dd = dd, dl = dl,
       sa = sa, sb = sb)
}

test_that("the parts and their standard errors follow the formulas", {
  r <- oaxbd(od$sen, od[c("sc", "pr")], od$g, ref = "B")
  h <- byhand(od$sen, od[c("sc", "pr")], od$g, "B")
  expect_identical(r$groups, c(a = "A", b = "B"))
  e <- setNames(r$dec$est, r$dec$part)
  s <- setNames(r$dec$se, r$dec$part)
  for (k in c("gap", "exp", "unexp", "expa", "unexpa", "inter")) {
    expect_equal(e[[k]], h[[k]], tolerance = 1e-12, info = k)
  }
  for (k in names(h$se)) expect_equal(s[[k]], h$se[[k]], tolerance = 1e-10)
  expect_equal(s[["gap"]], h$segap, tolerance = 1e-12)
  expect_identical(e[["endow"]], e[["exp"]])
  expect_identical(e[["coef"]], e[["unexpa"]])
  expect_equal(r$dec$lo, r$dec$est - qnorm(0.975) * r$dec$se,
               tolerance = 1e-15)
  expect_equal(r$coef$a, unname(h$ba), tolerance = 1e-12)
  expect_equal(r$coef$bse, unname(sqrt(diag(h$vb))), tolerance = 1e-10)
  expect_equal(r$means$b, unname(h$xb), tolerance = 1e-14)
  expect_identical(r$coef$term, c("const", "sc", "pr"))
  # by covariate
  b <- r$by
  expect_equal(b$exp, unname(h$dd * h$bb), tolerance = 1e-12)
  expect_equal(b$unexp, unname(h$xa * h$dl), tolerance = 1e-12)
  expect_equal(b$shift, unname(h$dl), tolerance = 1e-12)
  expect_equal(b$expse, unname(sqrt(h$dd^2 * diag(h$vb) +
                                      h$bb^2 * diag(h$sa + h$sb))),
               tolerance = 1e-10)
  expect_equal(b$unexpse, unname(sqrt(h$xa^2 * diag(h$va + h$vb) +
                                        h$dl^2 * diag(h$sa))),
               tolerance = 1e-10)
  expect_identical(b$expse[1], 0)
})

test_that("the identities hold exactly, and the choice of reference", {
  set.seed(8)
  n <- 300
  g <- sample(c("p", "q"), n, TRUE)
  x <- cbind(u = rnorm(n, g == "p"), v = runif(n))
  y <- 1 + x %*% c(0.5, -1) + 0.7 * (g == "p") * x[, 1] + rnorm(n)
  r <- oaxbd(y[, 1], x, g, ref = "q")
  e <- setNames(r$dec$est, r$dec$part)
  expect_equal(e[["gap"]], e[["exp"]] + e[["unexp"]], tolerance = 1e-13)
  expect_equal(e[["gap"]], e[["expa"]] + e[["unexpa"]], tolerance = 1e-13)
  expect_equal(e[["gap"]], e[["endow"]] + e[["coef"]] + e[["inter"]],
               tolerance = 1e-13)
  # the two explained parts differ by the interaction
  expect_equal(e[["expa"]] - e[["exp"]], e[["inter"]], tolerance = 1e-13)
  expect_equal(sum(r$by$exp), e[["exp"]], tolerance = 1e-13)
  expect_equal(sum(r$by$unexp), e[["unexp"]], tolerance = 1e-13)
  # the other reference: the gap changes sign, the parts swap roles
  s <- oaxbd(y[, 1], x, g, ref = "p")
  f <- setNames(s$dec$est, s$dec$part)
  expect_identical(s$groups, c(a = "q", b = "p"))
  expect_equal(f[["gap"]], -e[["gap"]], tolerance = 1e-14)
  expect_equal(f[["exp"]], -e[["expa"]], tolerance = 1e-13)
  expect_equal(f[["unexp"]], -e[["unexpa"]], tolerance = 1e-13)
  expect_identical(oaxbd(y[, 1], x, g, ref = 2), r)
  # moving a covariate's origin: the unexplained total stays, and c times
  # its shift moves between its line and the constant's
  x2 <- x
  x2[, "u"] <- x2[, "u"] + 3
  r2 <- oaxbd(y[, 1], x2, g, ref = "q")
  expect_equal(r2$dec$est, r$dec$est, tolerance = 1e-12)
  expect_equal(r2$by$unexp[2] - r$by$unexp[2], 3 * r$by$shift[2],
               tolerance = 1e-12)
  expect_equal(r2$by$unexp[1] - r$by$unexp[1], -3 * r$by$shift[2],
               tolerance = 1e-11)
})

test_that("the values agree with the earlier package's", {
  r <- oaxbd(od$sen, od[c("sc", "pr")], od$g, ref = "B")
  e <- setNames(r$dec$est, r$dec$part)
  o <- ref$oax$o
  expect_equal(e[["gap"]], o$gap, tolerance = 1e-13)
  expect_equal(unname(e[c("exp", "unexp")]), unname(o$twofold_B),
               tolerance = 1e-12)
  expect_equal(unname(e[c("expa", "unexpa")]), unname(o$twofold_A),
               tolerance = 1e-12)
  expect_equal(unname(e[c("endow", "coef", "inter")]), unname(o$threefold),
               tolerance = 1e-12)
  expect_equal(r$by$exp, o$by_variable$explained_B, tolerance = 1e-12)
  expect_equal(r$by$unexp, o$by_variable$unexplained_A, tolerance = 1e-12)
  expect_equal(r$by$shift, o$by_variable$unexplained_shift, tolerance = 1e-12)
})

test_that("edge cases: no covariates, one, and exact fits", {
  y <- c(1, 2, 3, 6, 7, 8)
  g <- c(1, 1, 1, 2, 2, 2)
  r <- oaxbd(y, NULL, g, ref = 2)
  # with the constant only, the whole gap is unexplained
  expect_equal(r$dec$est[1:3], c(-5, 0, -5), tolerance = 1e-14)
  expect_identical(r$by$term, "const")
  r <- oaxbd(y, c(1, 2, 3, 2, 3, 4), g)
  expect_identical(r$coef$term, c("const", "x"))
  # an exact fit in each group: no residual variance, standard errors of
  # the coefficients 0
  expect_equal(r$coef$ase, c(0, 0), tolerance = 1e-12)
  expect_identical(oaxbd(y, matrix(c(1, 2, 3, 2, 3, 4)), g)$coef$term,
                   c("const", "x"))
  y8 <- c(1, 2, 3, 5, 6, 7, 8, 10)
  x8 <- list(1:8, c(0, 1, 0, 1, 1, 0, 1, 0))
  expect_identical(oaxbd(y8, x8, rep(1:2, each = 4))$coef$term,
                   c("const", "x1", "x2"))
})

test_that("oaxbd() refuses what it cannot use", {
  y <- c(1, 2, 3, 4, 6, 7, 8, 9)
  g <- rep(1:2, each = 4)
  errs("oaxbd", list(
    list("exactly two groups", y, NULL, rep(1:4, 2)),
    list("code a factor", y, factor(rep(c("a", "b"), 4)), g),
    list("code a factor", y, data.frame(f = letters[1:8]), g),
    list("code a factor", y, 1:7, g),
    list("code a factor", y, matrix(1:14, 7), g),
    list("code a factor", y, "a", g),
    list("finite numbers", y, c(1:7, NA), g),
    list("more cases than coefficients, 4", y,
         cbind(1:8, c(1, 4, 2, 5, 3, 6, 4, 8), c(2, 1, 1, 3, 2, 2, 5, 1)), g),
    list("constant, or a combination", y, rep(1, 8), g),
    list("constant, or a combination", y, cbind(1:8, 2 * (1:8)), g),
    list("position from 1 to 2", y, NULL, g, ref = 3)
  ))
})
