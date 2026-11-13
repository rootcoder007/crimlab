# Cross-type K, the segregation test, quadrats, Morisita's index and the
# segregation indices. The kernels are checked in src/t_ppx.cpp; here the
# binding against spatstat and vegan (stored, versions named) and the
# formulas written out.

pdat <- function() {
  set.seed(101)
  n <- 80
  p <- cbind(x = runif(n, 0, 2), y = runif(n, 0, 1))
  m <- sample(c("a", "b", "c"), n, TRUE)
  list(p = p, m = m)
}

test_that("cross-type K and the segregation statistic agree with spatstat", {
  d <- pdat()
  k <- crossk(d$p, marks = d$m, i = "a", j = "b", win = c(0, 2, 0, 1),
              r = c(0, 0.05, 0.1, 0.2, 0.3), sigma = 0.15, nsim = 19, seed = 1)
  expect_s3_class(k, "crossk")
  expect_named(k, c("k", "seg", "i", "j"))
  # spatstat.explore 3.8.3: Kcross(X, "a", "b", r = , correction =
  #   "isotropic")$iso
  expect_equal(k$k$k, c(0, 0.00625, 0.0316896943210901, 0.116373975658453,
                        0.282281120265239), tolerance = 1e-13)
  expect_equal(k$k$l, sqrt(k$k$k / pi), tolerance = 1e-15)
  expect_equal(k$k$theo, pi * k$k$r^2, tolerance = 1e-15)
  # the same, with type c about itself
  kc <- crossk(d$p, marks = factor(d$m), i = 3, j = "c", win = c(0, 2, 0, 1),
               r = c(0, 0.1, 0.25))
  expect_equal(kc$k$k, c(0, 0.0174518997874244, 0.181548788198721),
               tolerance = 1e-13)
  expect_null(kc$seg)
  # the statistic of segregation.test(), sigma 0.15
  expect_equal(k$seg[["t"]], 6.50639406923695, tolerance = 1e-13)
  expect_identical(k, crossk(d$p, marks = d$m, i = "a", j = "b",
                             win = c(0, 2, 0, 1), r = c(0, 0.05, 0.1, 0.2, 0.3),
                             sigma = 0.15, nsim = 19, seed = 1))
  expect_gte(k$seg[["p"]], 0.05)
  expect_lte(k$seg[["p"]], 1)
  expect_output(print(k), "segregation T")
  expect_output(print(kc), "crossk: K of c about c")
  expect_true(is.nan(crossk(d$p, marks = d$m, i = "a", j = "b",
                            win = c(0, 2, 0, 1), r = 0.1, sigma = 0.15,
                            nsim = 0)$seg[["p"]]))
})

test_that("quadrat counts and Morisita's index agree with spatstat and vegan", {
  d <- pdat()
  q <- qdrtst(d$p, win = c(0, 2, 0, 1), nx = 4, ny = 2)
  expect_s3_class(q, "qdrtst")
  # spatstat.explore 3.8.3: quadrat.test(X, nx = 4, ny = 2) (two-sided)
  # (its counts as a vector run along the rows from the top)
  expect_identical(as.vector(t(q$counts)), c(9, 6, 12, 10, 9, 11, 10, 13))
  expect_equal(unname(q$stat[c("x2", "df", "ptwo")]),
               c(3.2, 7, 0.268190516527803), tolerance = 1e-13)
  expect_equal(q$stat[["p"]], pchisq(3.2, 7, lower.tail = FALSE),
               tolerance = 1e-15)
  expect_identical(q$stat[["lambda"]], 40)
  expect_output(print(q), "qdrtst: X2 3.2 on 7 df")
  # vegan 2.7.5: dispindmorisita(matrix(cnt, ncol = 1))
  m <- morsta(c(3, 0, 7, 12, 1, 4, 0, 9))
  expect_equal(unname(m[c("imor", "mclu", "muni", "imst", "p")]),
               c(1.67619047619048, 1.25750755070369, 0.848281976590782,
                 0.531048082636746, 7.1619349900267e-05), tolerance = 1e-13)
  e <- morsta(c(4, 5, 4, 6, 5, 4))
  expect_equal(unname(e[c("imor", "mclu", "muni", "imst", "p")]),
               c(0.841269841269841, 1.29009266644556, 0.845600430129136,
                 -0.502560659092045, 0.982175450874206), tolerance = 1e-13)
  # the other two arms of the standardised index
  s <- c(0, 0, 0, 0, 30)
  ms <- morsta(s)
  expect_equal(ms[["imst"]], 0.5 + 0.5 * (ms[["imor"]] - ms[["mclu"]]) /
                 (5 - ms[["mclu"]]), tolerance = 1e-15)
  mi <- morsta(c(0, 2, 5, 1, 3, 2))
  expect_equal(mi[["imst"]], 0.5 * (mi[["imor"]] - 1) / (mi[["mclu"]] - 1),
               tolerance = 1e-15)
  mu <- morsta(c(7, 13, 9, 11, 10, 12, 8, 10))
  expect_equal(mu[["imst"]], -0.5 * (mu[["imor"]] - 1) / (mu[["muni"]] - 1),
               tolerance = 1e-15)
  expect_identical(morsta(q$counts), morsta(as.vector(q$counts)))
})

test_that("segregation indices follow their formulas", {
  x <- cbind(a = c(90, 80, 10, 5, 50, 0), b = c(10, 20, 85, 95, 50, 0),
             c = c(5, 0, 30, 2, 10, 0))
  s <- dissix(x, delta = 0.3)
  expect_s3_class(s, "dissix")
  tt <- colSums(x)
  dd <- function(a, b) sum(abs(x[, a] / tt[a] - x[, b] / tt[b])) / 2
  expect_equal(s$d["a", "b"], dd("a", "b"), tolerance = 1e-15)
  expect_equal(s$d["b", "c"], dd("b", "c"), tolerance = 1e-15)
  expect_identical(diag(s$d), c(a = 0, b = 0, c = 0))
  y <- x[rowSums(x) > 0, ]
  t <- rowSums(y)
  tot <- sum(y)
  ent <- function(p) {
    ifelse(p > 0, -p * log(p), 0) + ifelse(p < 1, -(1 - p) * log(1 - p), 0)
  }
  for (k in 1:3) {
    P <- sum(y[, k]) / tot
    p <- y[, k] / t
    o <- (t - y[, k]) / (tot - sum(y[, k]))
    expect_equal(s$grp$is[k], sum(abs(y[, k] / sum(y[, k]) - o)) / 2,
                 tolerance = 1e-14)
    expect_equal(s$grp$h[k], sum(t * (ent(P) - ent(p))) / (ent(P) * tot),
                 tolerance = 1e-14)
    expect_equal(s$grp$atk[k], 1 - P / (1 - P) *
                   sum((1 - p)^0.7 * p^0.3 * t / (P * tot))^(1 / 0.7),
                 tolerance = 1e-14)
    expect_equal(s$grp$gini[k], sum(outer(t, t) * abs(outer(p, p, "-"))) /
                   (2 * tot^2 * P * (1 - P)), tolerance = 1e-14)
  }
  expect_null(s$dm)
  w <- spwt(nb = list(2, c(1, 3), c(2, 4), c(3, 5), 4, 0), type = "contig",
            style = "raw")
  m <- dissix(data.frame(x), w = w)
  cm <- matrix(0, 5, 5)
  cm[cbind(w$i, w$j)] <- 1
  z <- y[, "a"] / (y[, "a"] + y[, "b"])
  expect_equal(m$dm["a", "b"], dd("a", "b") -
                 sum(cm * abs(outer(z, z, "-"))) / sum(cm), tolerance = 1e-14)
  expect_output(print(m), "Morrill")
  expect_output(print(s), "dissix: dissimilarity")
  expect_identical(colnames(dissix(unname(x))$d), c("g1", "g2", "g3"))
})

test_that("the point-pattern functions refuse bad input in words", {
  d <- pdat()
  errs("crossk", list(
    list("at least 2 points", d$p[1, , drop = FALSE], marks = "a", i = "a",
         j = "a", win = c(0, 2, 0, 1), r = 0.1),
    list("inside `win`", d$p, marks = d$m, i = "a", j = "b",
         win = c(0, 1, 0, 1), r = 0.1),
    list("`r` must increase", d$p, marks = d$m, i = "a", j = "b",
         win = c(0, 2, 0, 1), r = c(0.2, 0.1)),
    list("names no type", d$p, marks = d$m, i = "z", j = "b",
         win = c(0, 2, 0, 1), r = 0.1),
    list("type `i` needs a point", d$p[1:2, ], marks = c("a", "b"), i = "a",
         j = "a", win = c(0, 2, 0, 1), r = 0.1),
    list("widen `sigma`", d$p, marks = d$m, i = "a", j = "b",
         win = c(0, 2, 0, 1), r = 0.1, sigma = 1e-4)
  ))
  errs("qdrtst", list(
    list("inside `win`", d$p, win = c(0, 1, 0, 1), nx = 2),
    list("from 2 to 1e7", d$p, win = c(0, 2, 0, 1), nx = 1)
  ))
  errs("morsta", list(
    list("at least 2 quadrats", 3),
    list("more than 1", c(1, 0)),
    list("`counts` must hold finite numbers from 0", c(1, -1)),
    list("`counts` must be numeric", c("a", "b")),
    list("`crit`", 1:3, crit = 1)
  ))
  errs("dissix", list(
    list("a matrix or data frame", 1:3),
    list("at least 2 groups", matrix(1:3)),
    list("positive total", cbind(a = 1:3, b = 0)),
    list("`delta`", cbind(1:3, 3:1), delta = 1)
  ))
})
