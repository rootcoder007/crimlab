# Prospective hot-spotting. The kernel is checked by hand in
# src/t_prh.cpp; here the binding against the sum written out in R.

test_that("the risk surface and its evaluation follow the formula", {
  set.seed(6)
  n <- 150
  p <- cbind(x = runif(n, 0, 1000), y = runif(n, 0, 800))
  day <- sample(1:120, n, TRUE)
  r <- prohs(
    p,
    t = day,
    at = 100,
    win = c(0, 1000, 0, 800),
    cell = 100,
    bw = 300,
    tw = 30,
    ahead = 14,
    top = 0.2
  )
  expect_s3_class(r, "prohs")
  expect_equal(r$x, seq(50, 950, 100))
  expect_equal(r$y, seq(50, 750, 100))
  use <- day < 100 & day >= 70
  z <- outer(
    r$x,
    r$y,
    Vectorize(function(cx, cy) {
      d <- sqrt((cx - p[use, 1])^2 + (cy - p[use, 2])^2)
      sum((d <= 300) / (1 + d / 100) / (1 + (100 - day[use]) / 7))
    })
  )
  expect_equal(r$z, z, tolerance = 1e-13)
  k <- ceiling(0.2 * length(z))
  expect_identical(sum(r$hot), as.integer(k))
  expect_true(min(r$z[r$hot]) >= max(r$z[!r$hot]))
  fut <- day >= 100 & day < 114
  ci <- floor(p[fut, 1] / 100) + 1
  cj <- floor(p[fut, 2] / 100) + 1
  h <- r$hot[cbind(ci, cj)]
  expect_equal(r$eval[["n"]], sum(fut))
  expect_equal(r$eval[["hits"]], sum(h))
  expect_equal(r$eval[["pai"]], mean(h) / (k / length(z)), tolerance = 1e-14)
  expect_identical(r$par[["bw"]], 300)
  expect_output(print(r), "10 x 8 grid")
  expect_identical(
    prohs(p, t = day, at = 100, win = c(0, 1000, 0, 800), cell = 100)$par[[
      "bw"
    ]],
    800
  )
})

test_that("prohs() refuses bad input in words", {
  p <- cbind(x = c(1, 2), y = c(1, 2))
  errs(
    "prohs",
    list(
      list("`t`", p, t = 1, at = 2, win = c(0, 5, 0, 5), cell = 1),
      list("`cell`", p, t = c(1, 2), at = 2, win = c(0, 5, 0, 5), cell = 0),
      list(
        "`top` must be above 0",
        p,
        t = c(1, 2),
        at = 2,
        win = c(0, 5, 0, 5),
        cell = 1,
        top = 0
      ),
      list(
        "more than 1e8 cells",
        p,
        t = c(1, 2),
        at = 2,
        win = c(0, 1e6, 0, 1e6),
        cell = 1,
        bw = 2
      )
    )
  )
})
