# Kernel intensity surfaces. The kernel is checked in src/t_hot.cpp (the
# direct double sum, the integral, the corrections by hand); here the
# binding against MASS and spatstat, stored.

test_that("the grid surface is kde2d() times n, with its bandwidth", {
  set.seed(61)
  p <- cbind(x = runif(30, 0, 4), y = runif(30, 0, 3))
  k <- kdsurf(p, win = c(0, 4, 0, 3), nx = 5, ny = 4)
  expect_s3_class(k, "kdsurf")
  expect_named(k, c("x", "y", "z", "h", "regions"))
  # MASS 7.3.65: MASS::kde2d(p[, 1], p[, 2], n = c(5, 4),
  #   lims = c(0, 4, 0, 3))$z
  ref <- c(0.0285023824711698, 0.0696718646083049, 0.0434175381865766,
           0.012376936127817, 0.00362754972840817, 0.029669084776741,
           0.0513319463396978, 0.0678721139611184, 0.0725909359905949,
           0.0287352453678222, 0.0830005912709508, 0.120592616231417,
           0.088129048619979, 0.0519328160337222, 0.0226172423990241,
           0.054906301269498, 0.0494240477537658, 0.0252650882345638,
           0.0226652212694971, 0.0259927766737544)
  expect_equal(as.vector(k$z) / 30, ref, tolerance = 1e-13)
  expect_equal(dim(k$z), c(5L, 4L))
  expect_equal(unname(k$h), c(MASS::bandwidth.nrd(p[, 1]),
                              MASS::bandwidth.nrd(p[, 2])) / 4,
               tolerance = 1e-14)
  # regions written out: the fewest cells, densest first, holding the share
  z <- sort(as.vector(k$z), decreasing = TRUE)
  m <- vapply(c(0.5, 0.9), function(s) which(cumsum(z) >= s * sum(z))[1], 1)
  expect_equal(k$regions$level, z[m], tolerance = 1e-15)
  expect_equal(k$regions$area, m * 1 * 1, tolerance = 1e-14)
  expect_equal(k$regions$share, m / 20, tolerance = 1e-15)
  expect_output(print(k), "kdsurf: 5 x 4 grid")
  # weights, two bandwidths and the default extent, against the direct sum
  w <- runif(30)
  g <- kdsurf(p, w = w, h = c(0.3, 0.5), nx = 7, ny = 6, cover = 0.5)
  d <- outer(g$x, g$y, Vectorize(function(a, b) {
    sum(w * dnorm(a, p[, 1], 0.3) * dnorm(b, p[, 2], 0.5))
  }))
  expect_equal(g$z, d, tolerance = 1e-13)
  expect_equal(range(g$x), range(p[, 1]) + c(-0.9, 0.9), tolerance = 1e-15)
  expect_equal(kdsurf(p, h = 0.3, nx = 7, ny = 6)$z,
               kdsurf(p, h = c(0.3, 0.3), nx = 7, ny = 6)$z)
})

test_that("intensities at the points match density.ppp()'s corrections", {
  set.seed(61)
  p <- cbind(x = runif(30, 0, 4), y = runif(30, 0, 3))
  wn <- c(0, 4, 0, 3)
  # spatstat.explore 3.8.3: density.ppp(X, sigma = 0.5, at = "points",
  #   leaveoneout = TRUE, edge = , diggle = ), first six points
  ref <- list(
    none = c(1.93706104695838, 0.570549570224546, 2.90362941412435,
             1.90034871608384, 3.29399636414429, 3.46774851488923),
    uniform = c(3.03655243685482, 1.39668087541567, 2.95801630386269,
                5.66684325622259, 5.78471433434642, 4.29977528738161),
    diggle = c(3.14212442234532, 0.850510479949603, 3.16737195414136,
               2.8949806598689, 5.28052162335014, 5.27315666456572))
  for (e in names(ref)) {
    k <- kdsurf(p, h = 0.5, win = wn, edge = e, at = "points")
    expect_equal(k$z[1:6], ref[[e]], tolerance = 1e-13, info = e)
  }
  expect_null(k$regions)
  expect_output(print(k), "intensity at 30 points")
  a <- kdsurf(p, h = 0.5, at = cbind(x = c(1, 2), y = c(1, 1)))
  expect_equal(a$z, c(sum(dnorm(1, p[, 1], 0.5) * dnorm(1, p[, 2], 0.5)),
                      sum(dnorm(2, p[, 1], 0.5) * dnorm(1, p[, 2], 0.5))),
               tolerance = 1e-14)
  # a grid node with none of its kernel inside the window
  u <- kdsurf(p, h = 0.01, win = wn, nx = 3, ny = 3, edge = "uniform")
  expect_false(anyNA(u$z))
})

test_that("kdsurf() refuses bad input in words", {
  p <- cbind(x = c(1, 2, 3), y = c(1, 3, 2))
  errs("kdsurf", list(
    list("at least 2 points", p[1, , drop = FALSE]),
    list("spread along both axes", cbind(x = 1:3, y = 1)),
    list("one or two positive numbers", p, h = c(1, 1, 1)),
    list("needs `win`", p, h = 1, edge = "uniform"),
    list("`at` must be NULL", p, h = 1, at = "grid"),
    list("at most 1e8", p, h = 1, nx = 1e5, ny = 1e5),
    list("shares between 0 and 1", p, h = 1, cover = 1),
    list("none of its kernel", p, h = 0.1, win = c(10, 11, 10, 11),
         edge = "diggle", at = "points")
  ))
})
