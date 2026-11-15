# Internal benchmarking. The kernel against glm() weights and p.adjust().

test_that("benchmarks follow the propensity weights of glm()", {
  set.seed(5)
  n <- 600
  off <- sample(c("a", "b", "c", "d"), n, TRUE)
  night <- rbinom(n, 1, ifelse(off == "a", 0.7, 0.3))
  pl <- runif(n)
  y <- rbinom(n, 1, plogis(-1 + night + 0.5 * (off == "b")))
  x <- data.frame(night, pl)
  r <- ibench(y, off, x, minn = 10)
  expect_s3_class(r, "ibench")
  for (o in c("a", "b", "c", "d")) {
    d <- as.numeric(off == o)
    e <- fitted(glm(
      d ~ night + pl,
      family = binomial,
      control = glm.control(epsilon = 1e-14, maxit = 100)
    ))
    w <- (e / (1 - e))[d == 0]
    b <- sum(w * y[d == 0]) / sum(w)
    ess <- sum(w)^2 / sum(w^2)
    rt <- mean(y[d == 1])
    z <- (rt - b) / sqrt(rt * (1 - rt) / sum(d) + b * (1 - b) / ess)
    k <- r$tab$off == o
    expect_equal(r$tab$bench[k], b, tolerance = 1e-8)
    expect_equal(r$tab$ess[k], ess, tolerance = 1e-8)
    expect_equal(r$tab$z[k], z, tolerance = 1e-7)
  }
  expect_equal(r$tab$q, p.adjust(r$tab$p, "BH"), tolerance = 1e-12)
  expect_identical(r$tab$flag, r$tab$q < 0.05)
  expect_output(print(r), "ibench: 4 officers, 4 tested")
})

test_that("officers with few stops are not tested; no x is the plain share", {
  y <- c(1, 0, 1, 1, 0, 0, 1, 0, 1)
  off <- c("a", "a", "a", "a", "b", "b", "b", "b", "c")
  r <- ibench(y, off, minn = 2)
  expect_true(is.na(r$tab$z[3]))
  expect_equal(r$tab$bench[1], mean(y[off != "a"]), tolerance = 1e-12)
  expect_identical(r$nx, 0L)
})

test_that("the doubly robust benchmark is the weighted regression by hand", {
  set.seed(8)
  n <- 800
  off <- sample(c("a", "b", "c", "d"), n, TRUE)
  night <- rbinom(n, 1, ifelse(off == "a", 0.7, 0.3))
  pl <- runif(n)
  y <- rbinom(n, 1, plogis(-1 + night + pl + 0.6 * (off == "b")))
  x <- data.frame(night, pl)
  r <- ibench(y, off, x, minn = 10, dr = TRUE)
  expect_true(r$dr)
  for (o in c("a", "b", "c", "d")) {
    d <- as.numeric(off == o)
    e <- fitted(glm(
      d ~ night + pl,
      family = binomial,
      control = glm.control(epsilon = 1e-14, maxit = 100)
    ))
    w <- ifelse(d == 1, 1, e / (1 - e))
    # the outcome regression with the propensity weights
    g <- suppressWarnings(glm(
      y ~ d + night + pl,
      family = quasibinomial,
      weights = w,
      control = glm.control(epsilon = 1e-14, maxit = 100)
    ))
    Z <- cbind(1, d, night, pl)
    m0 <- plogis(drop(Z[d == 1, ] %*% (coef(g) * c(1, 0, 1, 1))))
    m1 <- plogis(drop(Z[d == 1, ] %*% coef(g)))
    k <- r$tab$off == o
    expect_equal(r$tab$bench[k], mean(m0), tolerance = 1e-8)
    expect_equal(mean(m1), r$tab$rate[k], tolerance = 1e-8)
    # the sandwich variance of mean(m1 - m0)
    mu <- fitted(g)
    H <- crossprod(Z * (w * mu * (1 - mu)), Z)
    M <- crossprod(Z * (w * (y - mu)), Z * (w * (y - mu)))
    Zo <- Z[d == 1, ]
    Z0 <- Zo
    Z0[, 2] <- 0
    gr <- colMeans(Zo * (m1 * (1 - m1)) - Z0 * (m0 * (1 - m0)))
    v <- drop(t(gr) %*% solve(H, M) %*% solve(H, gr))
    expect_equal(r$tab$z[k], (mean(m1) - mean(m0)) / sqrt(v), tolerance = 1e-6)
  }
  # without x the doubly robust benchmark is the others' plain share
  y2 <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1)
  o2 <- rep(c("a", "b"), each = 5)
  r2 <- ibench(y2, o2, minn = 2, dr = TRUE)
  expect_equal(r2$tab$bench[1], mean(y2[o2 == "b"]), tolerance = 1e-12)
})

test_that("the local false discovery rate follows locfdr", {
  # the z statistics of 150 officers and the local rates the package
  # locfdr gives them with its defaults (written by dev/ref/lfdr.R)
  ref <- refget("lfdr")
  set.seed(7)
  no <- 150
  off <- rep(sprintf("o%03d", seq_len(no)), each = 40)
  night <- rbinom(no * 40, 1, rep(runif(no, 0.2, 0.8), each = 40))
  bad <- off %in% sprintf("o%03d", 1:8)
  y <- rbinom(no * 40, 1, plogis(-1 + night + 1.2 * bad))
  r <- ibench(y, off, data.frame(night), minn = 10, local = TRUE, fdr = 0.2)
  expect_equal(r$tab$z, ref$z, tolerance = 1e-10)
  expect_equal(r$tab$lfdr, ref$fdr, tolerance = 1e-10)
  expect_equal(unname(r$null), ref$null, tolerance = 1e-10)
  expect_named(r$null, c("delta", "sigma", "p0"))
  expect_identical(r$tab$flag, r$tab$lfdr < 0.2)
  expect_null(ibench(y, off, data.frame(night), minn = 10)$null)
})

test_that("ibench() refuses bad input in words", {
  errs(
    "ibench",
    list(
      list("`y` must be 0 or 1", c(2, 0), c("a", "b")),
      list("`off` must hold at least two", c(1, 0), c("a", "a")),
      list("`fdr`", c(1, 0), c("a", "b"), fdr = 0),
      list("`minn`", c(1, 0), c("a", "b"), minn = 0),
      list("`dr`", c(1, 0), c("a", "b"), dr = NA),
      list(
        "local false discovery rate could not be estimated",
        c(1, 0, 1, 1),
        c("a", "a", "b", "b"),
        minn = 1,
        local = TRUE
      ),
      list(
        "`x` must be a numeric matrix",
        c(1, 0),
        c("a", "b"),
        x = matrix(1, 3, 1)
      )
    )
  )
})

test_that("covariates as a matrix or one vector fit as the data frame", {
  set.seed(5)
  n <- 600
  off <- sample(c("a", "b", "c", "d"), n, TRUE)
  night <- rbinom(n, 1, ifelse(off == "a", 0.7, 0.3))
  pl <- runif(n)
  y <- rbinom(n, 1, plogis(-1 + night + 0.5 * (off == "b")))
  expect_identical(unclass(ibench(y, off, cbind(night, pl), minn = 10)),
                   unclass(ibench(y, off, data.frame(night, pl), minn = 10)))
  # one vector is the one-column frame of it
  expect_equal(unclass(ibench(y, off, night, minn = 10)),
               unclass(ibench(y, off, data.frame(night), minn = 10)),
               tolerance = 0)
  expect_error(ibench(y, off, cbind(night, c(pl[-1], NA)), minn = 10),
               "ibench: `x` must hold finite numbers with no missing values.",
               fixed = TRUE)
})

test_that("more columns of x than the work arrays allow are refused", {
  # (p + 2)^2 cells for each regression: past 8190 columns refused before
  # any of them is made; the input itself is 8 x 8191 zeros
  x <- matrix(0, 8, 8191)
  expect_error(ibench(rep(0:1, 4), rep(c("a", "b"), each = 4), x),
               "ibench: `x` must have at most 8190 columns.", fixed = TRUE)
})
