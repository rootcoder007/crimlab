# Regression discontinuity (rdloc(), rddens(), R/rdloc.R). The reference
# values come from base R alone, by the routes below written out from the
# papers: the bandwidth of Imbens and Kalyanaraman (2012, steps 1-3) with
# var(), lm() and the triangle's C_K = 480^(1/5); the local fits by
# weighted least squares with solve() and the hc1 sandwich by hand (the
# fuzzy one as the just-identified weighted IV); McCrary (2008) with its
# bins counted and its local-linear fits by lm(weights = ). The data are
# simulated; no number comes from the papers. Tolerance 1e-8 (1e-6 for
# the bandwidth, whose pilot fits are separate least-squares solves).

simrd <- function() {
  set.seed(20261120)
  n <- 1500
  d <- data.frame(r = runif(n, -1, 1))
  d$y <- 1 + 0.8 * d$r + 2 * d$r^2 * (d$r >= 0) + 0.5 * (d$r >= 0) +
    rnorm(n, 0, 0.3)
  d$t <- as.numeric(runif(n) < 0.2 + 0.6 * (d$r >= 0))
  d$y2 <- 1 + 0.8 * d$r + 0.7 * d$t + rnorm(n, 0, 0.3)
  d$g <- sample(1:60, n, TRUE)
  d
}

# Imbens and Kalyanaraman (2012), steps 1-3, written out from the paper
ikref <- function(y, x, c, ck) {
  n <- length(x)
  h1 <- 1.84 * stats::sd(x) * n^(-1 / 5)
  l <- x >= c - h1 & x < c
  r <- x >= c & x <= c + h1
  f <- (sum(l) + sum(r)) / (2 * n * h1)
  vl <- stats::var(y[l])
  vr <- stats::var(y[r])
  t <- x - c
  g <- stats::coef(stats::lm(y ~ I(x >= c) + t + I(t^2) + I(t^3)))
  m3 <- 6 * g[5]
  h2r <- 3.56 * (vr / (f * m3^2))^(1 / 7) * sum(x >= c)^(-1 / 7)
  h2l <- 3.56 * (vl / (f * m3^2))^(1 / 7) * sum(x < c)^(-1 / 7)
  sr <- x >= c & x <= c + h2r
  sl <- x >= c - h2l & x < c
  m2r <- 2 * stats::coef(stats::lm(y[sr] ~ t[sr] + I(t[sr]^2)))[3]
  m2l <- 2 * stats::coef(stats::lm(y[sl] ~ t[sl] + I(t[sl]^2)))[3]
  rr <- 2160 * vr / (sum(sr) * h2r^4)
  rl <- 2160 * vl / (sum(sl) * h2l^4)
  den <- f * ((m2r - m2l)^2 + rr + rl)
  unname(ck * ((vl + vr) / den)^(1 / 5) * n^(-1 / 5))
}

# the local fits by weighted least squares, the hc1 sandwich by hand
fitref <- function(y, r, h, p, tt = NULL) {
  w <- pmax(1 - abs(r / h), 0)
  k <- w > 0
  dd <- as.numeric(r >= 0)
  pw <- sapply(1:p, function(j) r^j)
  x <- cbind(1, dd, pw, dd * pw)[k, , drop = FALSE]
  w <- w[k]
  z <- x
  if (!is.null(tt)) x[, 2] <- tt[k]
  a <- solve(crossprod(z, w * x))
  b <- a %*% crossprod(z, w * y[k])
  e <- drop(y[k] - x %*% b)
  v <- a %*% crossprod(z * (w * e)) %*% t(a) * sum(k) / (sum(k) - 1)
  c(b[2], sqrt(v[2, 2]))
}

# McCrary (2008), written out from the paper, with given b and h
mcref <- function(r, c, b, h) {
  g <- floor((r - c) / b) * b + b / 2 + c
  xj <- seq(min(g), max(g), by = b)
  yj <- sapply(xj, function(v) sum(abs(g - v) < b / 4)) / (length(r) * b)
  fe <- sapply(c(FALSE, TRUE), function(right) {
    s <- if (right) xj > c else xj < c
    u <- xj[s] - c
    stats::coef(stats::lm(yj[s] ~ u, weights = pmax(1 - abs(u / h), 0)))[1]
  })
  c(log(fe[2]) - log(fe[1]),
    sqrt(24 / 5 / (length(r) * h) * (1 / fe[2] + 1 / fe[1])))
}

test_that("rdloc: the IK bandwidth, the two fits and their errors", {
  d <- simrd()
  f <- rdloc(d, "y", "r")
  expect_equal(f$h, 0.357513114301047, tolerance = 1e-6)
  expect_equal(f$h, ikref(d$y, d$r, 0, 480^(1 / 5)), tolerance = 1e-8)
  expect_identical(f$coef$method, c("robust", "conventional"))
  ref1 <- fitref(d$y, d$r, f$h, 1)
  ref2 <- fitref(d$y, d$r, f$h, 2)
  expect_equal(f$coef$est, c(ref2[1], ref1[1]), tolerance = 1e-8)
  expect_equal(f$coef$se, c(ref2[2], ref1[2]), tolerance = 1e-8)
  expect_equal(f$coef$est[2], 0.444074781321117, tolerance = 1e-6)
  g <- rdloc(d, "y2", "r", trt = "t", bw = "user", h = 0.4,
             ci = "conventional")
  expect_equal(c(g$coef$est[1], g$coef$se[1]),
               c(0.735200363501564, 0.0744309194340123), tolerance = 1e-8)
  expect_output(print(f), "Imbens-Kalyanaraman")
  expect_identical(summary(f), f$coef)
})

test_that("rdloc: placebo cutoffs use one side only; errors", {
  d <- simrd()
  f <- rdloc(d, "y", "r", bw = "user", h = 0.3, plac = c(-0.5, 0.5))
  s <- d$r < 0
  p1 <- fitref(d$y[s], d$r[s] + 0.5, 0.3, 1)
  expect_equal(c(f$placebo$est[1], f$placebo$se[1]), p1, tolerance = 1e-8)
  fc <- rdloc(d, "y", "r", bw = "user", h = 0.3, se = "cluster", clus = "g")
  expect_identical(fc$ncl, length(unique(d$g[abs(d$r) < 0.3])))
  expect_error(rdloc(d, "y", "r", bw = "user"), "needs `h`")
  expect_error(rdloc(d, "y", "r", h = 0.3), "with `bw")
  expect_error(rdloc(d, "y", "r", cut = 5), "both sides")
  expect_error(rdloc(d, "y", "r", bw = "user", h = 0.001), "fewer than")
  expect_error(rdloc(d, "y", "r", plac = 0), "must not hold")
})

test_that("rddens is McCrary's two steps", {
  d <- simrd()
  m <- rddens(d$r, b = 0.05, h = 0.5)
  expect_equal(c(m$tab$theta, m$tab$se),
               c(-0.0474225013558036, 0.161390621033137), tolerance = 1e-8)
  expect_equal(c(m$tab$theta, m$tab$se), unname(mcref(d$r, 0, 0.05, 0.5)),
               tolerance = 1e-8)
  a <- rddens(d$r)
  expect_equal(a$b, 2 * stats::sd(d$r) / sqrt(nrow(d)), tolerance = 1e-12)
  expect_true(a$h > 0)
  expect_error(rddens(c(1, NA)), "rddens: ")
  expect_output(print(m), "McCrary")
})
