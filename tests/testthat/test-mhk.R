# etasf(), mvhwk(), sthwk(): each log-likelihood against a direct sum
# written here, the fits against stats::optim() on that sum, recovery of
# known parameters from hwksim(), and the refusals.

ll_etas <- function(t, m, hor, th) {
  lam <- vapply(seq_along(t), function(i) {
    j <- which(t < t[i])
    th[1] + th[2] * sum(exp(th[3] * m[j]) * (t[i] - t[j] + th[4])^(-th[5]))
  }, 0)
  q <- th[5] - 1
  sum(log(lam)) - th[1] * hor -
    th[2] * sum(exp(th[3] * m) * (th[4]^(-q) - (hor - t + th[4])^(-q)) / q)
}

ll_mv <- function(t, ty, hor, mu, a, b) {
  lam <- vapply(seq_along(t), function(i) {
    j <- which(t < t[i])
    mu[ty[i]] + b * sum(a[ty[i], ty[j]] * exp(-b * (t[i] - t[j])))
  }, 0)
  g <- 1 - exp(-b * (hor - t))
  sum(log(lam)) - sum(mu) * hor - sum(colSums(a)[ty] * g)
}

ll_st <- function(t, x, y, hor, area, th) {
  lam <- vapply(seq_along(t), function(i) {
    j <- which(t < t[i])
    r2 <- (x[i] - x[j])^2 + (y[i] - y[j])^2
    th[1] + th[2] * th[3] * sum(exp(-th[3] * (t[i] - t[j])) *
                                  exp(-r2 / (2 * th[4]^2)) / (2 * pi * th[4]^2))
  }, 0)
  sum(log(lam)) - th[1] * hor * area - th[2] * sum(1 - exp(-th[3] * (hor - t)))
}

lx <- hwksim(c(log(0.3), 0.5, 0.6, 1), 800, kernel = "lomax", seed = 5)

test_that("ETAS: the direct sum, optim()'s maximum, the Lomax process", {
  set.seed(4)
  t4 <- lx$t[lx$t < 400]
  m <- rexp(length(t4), 2)
  f <- etasf(t4, 400, mag = m)
  th <- unname(f$par)
  dm <- m - min(m)
  expect_equal(f$loglik, ll_etas(t4, dm, 400, th), tolerance = 1e-10)
  o <- optim(c(log(th[1:2]), th[3], log(th[4]), log(th[5] - 1)),
             function(u) {
               -ll_etas(t4, dm, 400, c(exp(u[1:2]), u[3], exp(u[4]),
                                         1 + exp(u[5])))
             },
             method = "BFGS", control = list(reltol = 1e-14, maxit = 500))
  expect_lte(-o$value, f$loglik + 1e-6)
  expect_true(f$afree)
  # without marks: the Lomax Hawkes fit of hwkfit()
  e <- etasf(lx$t, 800)
  h <- hwkfit(lx$t, 800, kernel = "lomax")
  expect_equal(e$loglik, h$loglik, tolerance = 1e-6)
  expect_equal(e$branch, h$par[["eta"]], tolerance = 1e-3)
  expect_true(is.na(e$se[["alpha"]]))
  expect_output(print(e), "alpha fixed at 0")
  expect_identical(etasf(lx$t, 800, start = e$par, maxit = 5)$n, nrow(lx))
  expect_identical(etasf(lx$t, 800, mag = rep(2, nrow(lx)))$afree, FALSE)
})

test_that("multivariate: the direct sum, one type, optim()'s maximum", {
  a <- hwksim(c(log(0.3), 0.3, 1), 200, seed = 1)
  b <- hwksim(c(log(0.2), 0.4, 1), 200, seed = 2)
  d <- rbind(data.frame(t = a$t, type = "a"), data.frame(t = b$t, type = "b"))
  d <- d[order(d$t), ]
  f <- mvhwk(d$t, d$type, 200)
  ty <- match(d$type, c("a", "b"))
  expect_equal(f$loglik, ll_mv(d$t, ty, 200, f$mu, f$a, f$beta),
               tolerance = 1e-10)
  o <- optim(log(f$par), function(u) {
    p <- exp(u)
    -ll_mv(d$t, ty, 200, p[1:2], matrix(p[3:6], 2), p[7])
  }, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lte(-o$value, f$loglik + 1e-6)
  expect_named(f$par, c("mu.a", "mu.b", "a.a.a", "a.b.a", "a.a.b", "a.b.b",
                        "beta"))
  expect_equal(f$rad, max(Mod(eigen(f$a)$values)), tolerance = 1e-10)
  expect_output(print(f), "spectral radius")
  # one type is hwkfit()'s exponential process
  g <- mvhwk(a$t, rep("x", nrow(a)), 200)
  h <- hwkfit(a$t, 200)
  expect_equal(g$loglik, h$loglik, tolerance = 1e-7)
  expect_equal(g$beta, h$par[["beta"]], tolerance = 1e-3)
  fb <- mvhwk(d$t, d$type, 200, beta = 1, start = f$par[1:6])
  expect_identical(fb$beta, 1)
  expect_length(fb$par, 6L)
})

test_that("space-time: the direct sum, recovery, the intensity", {
  s <- hwksim(c(log(0.02), 0.5, 1), 300, win = c(0, 10, 0, 10), sigma = 0.3,
              seed = 3)
  f <- sthwk(s$t, s$x, s$y, 300, win = c(0, 10, 0, 10))
  expect_equal(f$loglik, ll_st(s$t, s$x, s$y, 300, 100, unname(f$par)),
               tolerance = 1e-10)
  z <- (f$par - c(0.02, 0.5, 1, 0.3)) / f$se
  expect_true(all(abs(z) < 4))
  q <- data.frame(t = c(50, 120), x = c(5, 2), y = c(5, 8))
  e <- sthwk(s$t, s$x, s$y, 300, win = c(0, 10, 0, 10), par = f$par, at = q)
  expect_false(e$fit)
  expect_equal(e$loglik, f$loglik, tolerance = 1e-12)
  th <- unname(f$par)
  want <- vapply(1:2, function(i) {
    j <- which(s$t < q$t[i])
    r2 <- (q$x[i] - s$x[j])^2 + (q$y[i] - s$y[j])^2
    th[1] + th[2] * th[3] * sum(exp(-th[3] * (q$t[i] - s$t[j])) *
                                  exp(-r2 / (2 * th[4]^2)) / (2 * pi * th[4]^2))
  }, 0)
  expect_equal(e$lam, want, tolerance = 1e-12)
  expect_output(print(e), "evaluated at the given parameters")
  bb <- sthwk(s$t, s$x, s$y, 300, start = f$par, maxit = 3)
  expect_equal(bb$area, diff(range(s$x)) * diff(range(s$y)))
  expect_named(hwksim(f, 50, seed = 1), c("t", "x", "y", "gen"))
})

test_that("bad input is refused", {
  errs("etasf", list(
    list("at least 10", 1:3, 5),
    list("`mag` must have length", lx$t, 800, mag = 1:3),
    list("p above 1", lx$t, 800, start = c(1, 1, 0, 1, 0.5))
  ))
  errs("mvhwk", list(
    list("at most 10", lx$t, seq_len(nrow(lx)) %% 11, 800),
    list("positive", lx$t, rep(1:2, length.out = nrow(lx)), 800,
         start = rep(0, 7))
  ))
  t5 <- 1:12
  errs("sthwk", list(
    list("on a line", t5, rep(1, 12), 1:12, 20),
    list("lie in `win`", t5, 1:12, 1:12, 20, win = c(0, 5, 0, 5)),
    list("columns t, x and y", t5, 1:12, 12:1, 20, at = 3),
    list("columns t, x and y", t5, 1:12, 12:1, 20, at = list(t = 1, x = 2)),
    list("one length", t5, 1:12, 12:1, 20, at = list(t = 1, x = 2:3, y = 1)),
    list("eta below 1", t5, 1:12, 12:1, 20, start = c(1, 2, 1, 1)),
    list("`par` must", t5, 1:12, 12:1, 20, par = c(1, 2, 1))
  ))
})
