# sdbsim() and sdbthr(): the stability analysis against eigen() of the
# linearised model written out here, the continuum scheme's one-step
# growth of a small wave against that matrix, the agent model's burglary
# count against its binomial mean, and the refusals.

# the linearisation of the continuum model at wavenumber-squared k
jac <- function(k, a0, om, eps, dif, gam, eta, z) {
  b <- eps * dif * gam / om
  a <- a0 + b
  r <- gam / a
  matrix(c(-om + eps * dif * r - eta * dif / z * k,
           -r + 2 * dif * r / (z * a) * k,
           eps * dif * a,
           -a - dif / z * k), 2)
}

test_that("the stability analysis is eigen()'s", {
  k <- c(0, 0.5, 2, 10, 48.5)
  s <- sdbthr(1, 1, 2, 1, 1, 0.01, k = k)
  for (i in seq_along(k)) {
    m <- jac(k[i], 1, 1, 2, 1, 1, 0.01, 1)
    expect_equal(s$disp$trace[i], sum(diag(m)), tolerance = 1e-14)
    expect_equal(s$disp$det[i], det(m), tolerance = 1e-12)
    expect_equal(s$disp$growth[i], max(Re(eigen(m)$values)),
                 tolerance = 1e-12)
  }
  expect_true(s$unstable)
  expect_equal(s$detcrit, det(jac(s$kcrit, 1, 1, 2, 1, 1, 0.01, 1)),
               tolerance = 1e-12)
  expect_equal(unname(s$steady), c(2, 3, 1 / 3),
               tolerance = testthat_tolerance())
  expect_output(print(s), "hot spots form")
  st <- sdbthr(1, 1, 0, 1, 1, 0.01)
  expect_false(st$unstable)
  expect_null(st$disp)
  expect_output(print(st), "stable")
  # z = 4 with the lattice's neighbour count on both diffusions
  s4 <- sdbthr(1, 1, 2, 1, 1, 0.01, z = 4, k = 1)
  expect_equal(s4$disp$det, det(jac(1, 1, 1, 2, 1, 1, 0.01, 4)),
               tolerance = 1e-12)
})

test_that("the continuum scheme grows a small wave as the matrix says", {
  n <- 16
  a0 <- 0.5
  om <- 1
  eps <- 2
  dif <- 1
  gam <- 0.3
  eta <- 0.2
  z <- 4
  dt <- 0.01
  bb <- eps * dif * gam / om
  rb <- gam / (a0 + bb)
  wv <- outer(cos(2 * pi * 2 * (0:(n - 1)) / n), rep(1, n))
  e <- c(1e-7, 3e-8)
  r <- sdbsim(n, 1, "pde", a0, om, eta, gam, dt = dt, b = bb + e[1] * wv,
              rho = rb + e[2] * wv, dif = dif, eps = eps, z = z)
  pr <- function(f, m) sum((f - m) * wv) / sum(wv^2)
  k <- 2 - 2 * cos(2 * pi * 2 / n)
  want <- (diag(2) + dt * jac(k, a0, om, eps, dif, gam, eta, z)) %*% e
  expect_equal(c(pr(r$b, bb), pr(r$rho, rb)), as.vector(want),
               tolerance = 1e-5)
  expect_output(print(r), "pde model")
})

test_that("the agent model: burglaries at their binomial mean, seeded", {
  rho <- matrix(200, 10, 10)
  a <- sdbsim(10, 1, a0 = 0.1, omega = 0.1, gam = 0, b = matrix(0.2, 10, 10),
              rho = rho, seed = 4)
  p <- 1 - exp(-0.3)
  m <- 20000 * p
  expect_lt(abs(a$burg - m), 4 * sqrt(20000 * p * (1 - p)))
  expect_identical(sum(a$tot), a$burg)
  expect_identical(a, sdbsim(10, 1, a0 = 0.1, omega = 0.1, gam = 0,
                             b = matrix(0.2, 10, 10), rho = rho, seed = 4))
  # no offenders: b falls by 1 - omega a step
  d <- sdbsim(5, 4, b = matrix(1, 5, 5), gam = 0, omega = 0.25, seed = 1)
  expect_equal(d$meanb, 0.75^(1:4), tolerance = testthat_tolerance())
  expect_output(print(d), "agent model")
  expect_output(print(sdbsim(5, 0, seed = 1)), "0 steps")
  # suppression
  s <- matrix(FALSE, 8, 8)
  s[3, 3] <- TRUE
  q <- sdbsim(8, 3, "pde", sup = s, rho = matrix(0.1, 8, 8))
  expect_identical(q$tot[3, 3], 0)
})

test_that("bad input is refused", {
  errs("sdbsim", list(
    list("`n`", 2, 1),
    list("at most 2000", 2001, 1),
    list("`mode` must be one of", 5, 1, mode = "lattice"),
    list("`eta`", 5, 1, eta = 2),
    list("at most 1", 5, 1, omega = 2),
    list("n x n matrix", 5, 1, b = matrix(0, 4, 4)),
    list("whole numbers", 5, 1, rho = matrix(0.5, 5, 5)),
    list("`sup` must be an n x n", 5, 1, sup = TRUE),
    list("`sup`", 5, 1, sup = matrix(1, 5, 5)),
    list("unstable", 5, 1, mode = "pde", dt = 2)
  ))
  errs("sdbthr", list(
    list("`a0`", 0, 1, 1, 1, 1, 1),
    list("`eps`", 1, 1, -1, 1, 1, 1),
    list("`k`", 1, 1, 1, 1, 1, 1, k = -1)
  ))
})
