# hypq(): the stationary law against the generator solved here with
# solve(), Erlang's loss and delay formulas, by-hand travel times; and
# sqrtlw() against its closed form and a simulation of nearest distances.

# the hypercube's generator, built and solved here
ctmc <- function(lam, pref, mu, wait = FALSE) {
  n <- ncol(pref)
  s <- 2^n
  busy <- function(b) bitwAnd(b, 2^(0:(n - 1))) > 0
  q <- matrix(0, s, s)
  rho <- sum(lam) / (n * mu)
  for (b in 0:(s - 1)) {
    bb <- busy(b)
    if (!all(bb)) {
      for (j in seq_along(lam)) {
        u <- pref[j, !bb[pref[j, ]]][1] - 1
        q[b + 1, b + 2^u + 1] <- q[b + 1, b + 2^u + 1] + lam[j]
      }
    }
    sc <- if (all(bb) && wait) 1 - rho else 1
    for (u in which(bb) - 1) {
      q[b + 1, b - 2^u + 1] <- q[b + 1, b - 2^u + 1] + mu * sc
    }
    q[b + 1, b + 1] <- -sum(q[b + 1, ])
  }
  a <- t(q)
  a[s, ] <- 1
  solve(a, c(rep(0, s - 1), 1))
}

pref <- rbind(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1), c(3, 1, 2))
lam <- c(0.4, 0.7, 0.3, 0.6)

test_that("the stationary law is the generator's, and Erlang's", {
  r <- hypq(lam, pref, mu = 1.5)
  p <- ctmc(lam, pref, 1.5)
  expect_equal(r$state$p, p, tolerance = 1e-12)
  expect_identical(r$state$busy[c(1, 2, 8)], c("000", "100", "111"))
  a <- sum(lam) / 1.5
  eb <- (a^3 / 6) / sum(a^(0:3) / factorial(0:3))
  expect_equal(r$pfull, eb, tolerance = 1e-12)
  expect_equal(sum(r$wl), a * (1 - eb), tolerance = 1e-12)
  expect_equal(colSums(r$disp), lam / sum(lam) * (1 - eb), tolerance = 1e-12)
  expect_identical(r$lq, 0)
  w <- hypq(lam, pref, mu = 1.5, queue = "wait")
  expect_equal(w$state$p, ctmc(lam, pref, 1.5, TRUE), tolerance = 1e-12)
  ec <- 3 * eb / (3 - a * (1 - eb))
  expect_equal(w$pfull, ec, tolerance = 1e-12)
  expect_equal(w$lq, ec * a / (3 - a), tolerance = 1e-12)
  expect_equal(w$wq, w$lq / sum(lam))
  expect_output(print(w), "mean queue")
})

test_that("orders as a list; travel times", {
  r <- hypq(lam, pref, mu = 1.5)
  l <- hypq(lam, split(pref, row(pref)), mu = 1.5)
  expect_identical(l$state, r$state)
  tt <- matrix(c(2, 5, 8, 6, 3, 7, 9, 6, 2, 7, 8, 3), 3)
  v <- hypq(lam, pref, mu = 1.5, tt = tt)
  expect_equal(v$travel$mean, sum(r$disp * tt) / sum(r$disp))
  expect_equal(v$travel$atom, colSums(r$disp * tt) / colSums(r$disp))
  expect_output(print(v), "mean travel time")
})

test_that("the square-root law: the closed form and a simulation", {
  s <- sqrtlw(100, 4)
  expect_identical(s$dist, 2.5)
  expect_true(is.na(s$time))
  g <- sqrtlw(100, 8, busy = 0.5, metric = "grid", speed = 2)
  expect_equal(g$dist, sqrt(2 * pi) / 4 * 5)
  expect_identical(g$avail, 4)
  expect_equal(g$time, g$dist / 2)
  # a Poisson pattern of 400 units on a 20 x 20 torus: the mean distance
  # from random calls to the nearest unit is 1/2 sqrt(400 / 400)
  set.seed(5)
  ux <- runif(400, 0, 20)
  uy <- runif(400, 0, 20)
  d <- vapply(seq_len(2000), function(i) {
    cx <- runif(1, 0, 20)
    cy <- runif(1, 0, 20)
    dx <- pmin(abs(ux - cx), 20 - abs(ux - cx))
    dy <- pmin(abs(uy - cy), 20 - abs(uy - cy))
    min(sqrt(dx^2 + dy^2))
  }, 0)
  expect_equal(mean(d), sqrtlw(400, 400)$dist, tolerance = 0.05)
})

test_that("bad input is refused", {
  errs("hypq", list(
    list("`lam`", c(1, 0), pref[1:2, ]),
    list("one row per atom", lam, pref[1:3, ]),
    list("one order per atom", lam, list(1:3)),
    list("list all 3 units", lam[1:2], list(1:3, 1:2)),
    list("from 1 to 3", lam, pref + 1),
    list("each unit once", lam, rbind(pref[1:3, ], c(1, 1, 2))),
    list("rho below 1", lam, pref, mu = 0.5, queue = "wait"),
    list("at most 10 units", 1, matrix(1:11, 1)),
    list("one row per unit", lam, pref, tt = matrix(1, 4, 4)),
    list("`queue` must be one of", lam, pref, queue = "inf"),
    list("`tt` must hold finite numbers", lam, pref, tt = matrix(-1, 3, 4)),
    list("`pref` must hold numbers", lam, matrix("a", 4, 3))
  ))
  errs("sqrtlw", list(
    list("`area`", 0, 1),
    list("below 1", 1, 1, busy = 1),
    list("`metric` must be one of", 1, 1, metric = "manhattan"),
    list("`speed`", 1, 1, speed = -1)
  ))
})
