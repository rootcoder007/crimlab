# Hot-spot accuracy. The kernel is checked by hand in src/t_hsa.cpp; here
# the binding against order() and sums in base R.

test_that("PAI, PEI and RRI follow the hot cells", {
  set.seed(12)
  risk <- rgamma(400, 0.5)
  before <- rpois(400, risk)
  after <- rpois(400, risk)
  pr <- before + runif(400) / 10
  r <- hsacc(pr, after, top = c(0.05, 0.1), prev = before)
  expect_s3_class(r, "data.frame")
  for (k in 1:2) {
    m <- c(20, 40)[k]
    hot <- order(-pr)[1:m]
    n <- sum(after[hot])
    expect_identical(r$cells[k], m)
    expect_equal(r$hit[k], n / sum(after), tolerance = 1e-15)
    expect_equal(r$pai[k], (n / sum(after)) / (m / 400), tolerance = 1e-14)
    top <- sum(sort(after, decreasing = TRUE)[1:m])
    expect_identical(r$nmax[k], as.numeric(top))
    expect_equal(r$pei[k], n / r$nmax[k], tolerance = 1e-15)
    expect_equal(
      r$rri[k],
      (n / sum(before[hot])) / (sum(after) / sum(before)),
      tolerance = 1e-14
    )
  }
  h <- hsacc(seq_len(400) %in% order(-pr)[1:20], after)
  expect_identical(h$cells, 20)
  expect_equal(h$pai, r$pai[1], tolerance = 1e-15)
  expect_true(is.na(h$rri))
})

test_that("unequal areas: cells are taken until their area reaches the share", {
  r <- hsacc(c(3, 2, 1, 0), c(1, 1, 1, 1), area = c(1, 3, 2, 4), top = 0.3)
  expect_identical(r$cells, 2)
  expect_identical(r$area, 4)
  expect_equal(r$pai, (2 / 4) / (4 / 10), tolerance = 1e-15)
})

test_that("tied scores at the cut-off give one result in every order", {
  set.seed(3)
  lam <- rgamma(400, 0.6, 0.6)
  pred <- rpois(400, lam)
  obs <- rpois(400, lam)
  prev <- rpois(400, lam)
  ar <- runif(400, 0.5, 2)
  top <- c(0.05, 0.1, 0.2)
  for (tk in c("in", "out", "frac")) {
    a <- hsacc(pred, obs, area = ar, top = top, prev = prev, ties = tk)
    for (s in 1:3) {
      o <- sample(400)
      b <- hsacc(pred[o], obs[o], area = ar[o], top = top, prev = prev[o],
                 ties = tk)
      expect_equal(b, a, tolerance = 1e-12, info = tk)
    }
  }
  # four cells tie at the top; 20% of ten cells is two of them
  p <- c(5, 5, 5, 5, 4, 3, 2, 1, 0, 0)
  o <- c(3, 0, 5, 1, 0, 0, 2, 9, 0, 0)
  r <- function(tk) hsacc(p, o, top = 0.2, ties = tk)
  expect_identical(r("in")$cells, 4)
  expect_identical(r("in")$n, 9)
  expect_identical(r("out")$cells, 0)
  expect_identical(r("frac")$cells, 2)
  expect_identical(r("frac")$n, 4.5)
  expect_identical(r("frac")$area, 2)
  # untied scores: every choice gives the earlier cells
  q <- c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
  for (tk in c("out", "frac")) {
    expect_identical(hsacc(q, o, top = c(0.2, 0.3), ties = tk),
                     hsacc(q, o, top = c(0.2, 0.3)))
  }
  expect_identical(hsacc(q, o, top = 0.2)$cells, 2)
})

test_that("nmax does not depend on the order of cells of equal density", {
  # three cells of density 1 after the hot one; the larger taken first
  a <- hsacc(c(4, 3, 2, 1), c(0, 2, 1, 2), area = c(3, 2, 1, 2), top = 0.375)
  b <- hsacc(c(4, 3, 1, 2), c(0, 2, 2, 1), area = c(3, 2, 2, 1), top = 0.375)
  expect_identical(a, b)
  expect_identical(a$nmax, 4)
})

test_that("hsacc() refuses bad input in words", {
  errs(
    "hsacc",
    list(
      list("`obs` has no crimes", c(1, 2), c(0, 0)),
      list("`top` must be at most 1", c(1, 2), c(1, 0), top = 2),
      list(
        "`pred` must have one value per cell",
        c(TRUE, FALSE, TRUE),
        c(1, 0)
      ),
      list("`area` must be greater than 0", c(1, 2), c(1, 0), area = c(0, 1)),
      list("`ties` must be one of", c(1, 2), c(1, 0), ties = "all")
    )
  )
})
