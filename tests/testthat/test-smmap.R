sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))

test_that("smmap() shares one set of breaks over the periods", {
  p <- list(sq(0, 0), sq(1, 0), sq(2, 0))
  v <- cbind(`2023` = c(1, 5, 9), `2024` = c(2, NA, 12))
  r <- drw(smmap(p, v, brk = 2, meth = "equal"))
  br <- seq(1, 12, length.out = 3)
  expect_identical(r$breaks, br)
  expect_identical(r$areas$prd, rep(c("2023", "2024"), each = 3))
  expect_identical(r$areas$area, rep(1:3, 2))
  expect_identical(r$areas$v, as.vector(v))
  cls <- findInterval(as.vector(v), br, rightmost.closed = TRUE,
                      all.inside = TRUE)
  expect_identical(r$areas$cls, cls)
  expect_identical(r$areas$col, blpal("sequential", 2)[cls])
  expect_identical(r$bbox, c(0, 3, 0, 1))
  q <- drw(smmap(p, v, brk = 4))
  qb <- unique(stats::quantile(v, seq(0, 1, length.out = 5), na.rm = TRUE,
                               names = FALSE))
  expect_equal(q$breaks, qb, tolerance = 1e-15)
  d <- drw(smmap(p, as.data.frame(v)))
  expect_identical(d$areas$prd, rep(c("2023", "2024"), each = 3))
  u <- drw(smmap(p, unname(v)))
  expect_identical(u$areas$prd, rep(c("1", "2"), each = 3))
  k <- keeps(smmap(p, v, main = "By year"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("smmap() refuses what it cannot draw", {
  p <- list(sq(0, 0), sq(1, 0))
  msg <- paste("smmap: `v` must be a numeric matrix, one row per area and",
               "one column per period.")
  expect_error(smmap(p, 1:2), msg, fixed = TRUE)
  expect_error(smmap(p, matrix(1, 3, 2)), msg, fixed = TRUE)
  expect_error(smmap(p, matrix("a", 2, 2)), msg, fixed = TRUE)
  expect_error(smmap(p, data.frame(a = c("x", "y"))), msg, fixed = TRUE)
  expect_error(smmap(list(), matrix(1)),
               "smmap: `polys` must be a non-empty list of areas.",
               fixed = TRUE)
  expect_error(smmap(p, matrix(NA_real_, 2, 2)),
               "smmap: `v` has no finite value.", fixed = TRUE)
  expect_error(smmap(p, matrix(1:4, 2), brk = 0),
               "smmap: `brk` must be a whole number", fixed = TRUE)
  expect_error(smmap(p, matrix(1:4, 2), main = c("a", "b")),
               "smmap: `main` must be one string or NULL.", fixed = TRUE)
})
