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

test_that("ibench() refuses bad input in words", {
  errs(
    "ibench",
    list(
      list("`y` must be 0 or 1", c(2, 0), c("a", "b")),
      list("`off` must hold at least two", c(1, 0), c("a", "a")),
      list("`fdr`", c(1, 0), c("a", "b"), fdr = 0),
      list("`minn`", c(1, 0), c("a", "b"), minn = 0),
      list(
        "`x` must be a numeric matrix",
        c(1, 0),
        c("a", "b"),
        x = matrix(1, 3, 1)
      )
    )
  )
})
