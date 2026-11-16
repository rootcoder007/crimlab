sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))

test_that("chomap() classes the areas by quantile, equal or given breaks", {
  p <- list(a = sq(0, 0), b = sq(1, 0), c = sq(0, 1), d = sq(1, 1),
            e = rbind(sq(2, 0), NA, sq(2, 1)))
  v <- c(3, 10, 4, 7, NA)
  r <- drw(chomap(p, v, brk = 3))
  br <- unique(stats::quantile(v, seq(0, 1, length.out = 4), na.rm = TRUE,
                               names = FALSE))
  expect_equal(r$breaks, br, tolerance = 1e-15)
  cls <- findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE)
  expect_identical(r$areas$cls, cls)
  expect_identical(r$areas$area, names(p))
  expect_identical(r$areas$v, v)
  expect_identical(r$areas$col, blpal("sequential", length(br) - 1L)[cls])
  expect_identical(r$bbox, c(0, 3, 0, 2))
  e <- drw(chomap(p, v, brk = 2, meth = "equal"))
  expect_identical(e$breaks, seq(3, 10, length.out = 3))
  expect_identical(e$areas$cls, c(1L, 2L, 1L, 2L, NA))
  g <- drw(chomap(p, v, brk = c(0, 5, 10)))
  expect_identical(g$areas$cls, findInterval(v, c(0, 5, 10),
                                             rightmost.closed = TRUE,
                                             all.inside = TRUE))
  # repeated quantile breaks are merged: 1, 1, 1, 2 in two classes
  q <- drw(chomap(p[1:4], c(1, 1, 1, 2), brk = 2))
  expect_identical(q$breaks, c(1, 2))
  expect_identical(q$areas$cls, rep(1L, 4))
  one <- drw(chomap(p[1], 2))
  expect_identical(one$breaks, c(1.5, 2.5))
  expect_identical(one$bbox, c(0, 1, 0, 1))
  u <- drw(chomap(unname(p[1:2]), c(1, 2), brk = 2,
                  col = c("red", "blue")))
  expect_identical(u$areas$area, 1:2)
  expect_identical(u$areas$col, c("red", "blue"))
  k <- keeps(chomap(p, v, main = "Rates", pal = "grey"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("chomap() refuses what it cannot draw", {
  p <- list(sq(0, 0), sq(1, 0))
  expect_error(chomap(list(), 1),
               "chomap: `polys` must be a non-empty list of areas.",
               fixed = TRUE)
  expect_error(chomap(data.frame(x = 1:3, y = 1:3), 1),
               "chomap: `polys` must be a non-empty list", fixed = TRUE)
  msg <- paste("chomap: `polys[[2]]` must be a two-column numeric matrix of",
               "at least three points (rows of NA between rings).")
  expect_error(chomap(list(sq(0, 0), sq(0, 0)[1:2, ]), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(list(sq(0, 0), sq(0, 0) * Inf), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(list(sq(0, 0), rbind(sq(0, 0), c(1, NA))), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(p, 1), "chomap: `v` must be numbers, one per area.",
               fixed = TRUE)
  expect_error(chomap(p, c("a", "b")), "chomap: `v` must be numbers",
               fixed = TRUE)
  expect_error(chomap(p, c(NA, Inf)), "chomap: `v` has no finite value.",
               fixed = TRUE)
  expect_error(chomap(p, 1:2, meth = "jenks"),
               "chomap: `meth` must be \"quantile\" or \"equal\".",
               fixed = TRUE)
  msg <- paste("chomap: `brk` must be a whole number of classes from 1 to 9,",
               "or increasing breaks covering `v`.")
  expect_error(chomap(p, 1:2, brk = 10), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = 2.5), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = c(1.5, 3)), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = c(2, 1)), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, main = 1),
               "chomap: `main` must be one string or NULL.", fixed = TRUE)
  expect_error(chomap(p, 1:2, pal = "neon"),
               "chomap: the palette must be one of", fixed = TRUE)
  expect_error(chomap(p, 1:2, col = "nocolour"),
               "chomap: `col` must be colours.", fixed = TRUE)
})
