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
      list("`area` must be greater than 0", c(1, 2), c(1, 0), area = c(0, 1))
    )
  )
})
