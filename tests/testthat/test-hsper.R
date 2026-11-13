# Hot-spot persistence. The kernel is checked by hand in src/t_hot.cpp;
# here the binding against the rule written out in R.

hp <- function(p, prd, cell, top) {
  ix <- floor((p[, 1] - min(p[, 1])) / cell)
  iy <- floor((p[, 2] - min(p[, 2])) / cell)
  key <- paste(ix, iy)
  u <- unique(data.frame(ix, iy))
  u <- u[order(u$ix, u$iy), ]
  cid <- match(key, paste(u$ix, u$iy))
  pf <- factor(prd)
  cnt <- unclass(table(factor(cid, seq_len(nrow(u))), pf))
  k <- max(1, ceiling(top * nrow(u)))
  hot <- apply(cnt, 2, function(c) {
    cut <- sort(c, decreasing = TRUE)[k]
    h <- c > cut | (c == cut & c > 0)
    if (sum(h) > 2 * k) h <- c > cut
    h
  })
  list(cid = cid, cnt = cnt, hot = hot, k = k)
}

test_that("hot cells, shares and overlaps follow the rule", {
  set.seed(71)
  n <- 500
  p <- cbind(x = c(rnorm(200, 3, 0.4), runif(300, 0, 10)),
             y = c(rnorm(200, 3, 0.4), runif(300, 0, 10)))
  prd <- sample(c("2022", "2023", "2024"), n, replace = TRUE)
  h <- hsper(p, prd = prd, cell = 1, top = 0.1)
  r <- hp(p, prd, 1, 0.1)
  expect_s3_class(h, "hsper")
  expect_named(h, c("periods", "cells", "hot", "cell", "k"))
  expect_identical(h$cell, r$cid)
  expect_identical(h$k, as.integer(r$k))
  expect_identical(h$hot, unname(r$hot))
  expect_identical(h$periods$period, c("2022", "2023", "2024"))
  expect_identical(h$periods$n, as.numeric(colSums(r$cnt)))
  expect_equal(h$periods$hotshr,
               unname(colSums(r$cnt * r$hot) / colSums(r$cnt)),
               tolerance = 1e-15)
  half <- apply(r$cnt, 2, function(c) {
    which(cumsum(sort(c, decreasing = TRUE)) >= sum(c) / 2)[1]
  })
  expect_identical(h$periods$nhalf, as.numeric(half))
  jac <- sum(r$hot[, 1] & r$hot[, 2]) / sum(r$hot[, 1] | r$hot[, 2])
  expect_equal(h$periods$jacc, c(NA, jac, h$periods$jacc[3]),
               tolerance = 1e-15)
  expect_identical(h$cells$n, as.numeric(rowSums(r$cnt)))
  expect_identical(h$cells$nhot, as.numeric(rowSums(r$hot)))
  run <- apply(r$hot, 1, function(v) {
    x <- rle(v)
    max(0, x$lengths[x$values])
  })
  expect_identical(h$cells$run, as.numeric(run))
  expect_equal(h$cells$x[h$cell], min(p[, 1]) +
                 (floor((p[, 1] - min(p[, 1]))) + 0.5), tolerance = 1e-14)
  expect_output(print(h), "hsper: .* active cells")
})

test_that("ties, empty hot sets, factors and longitude and latitude", {
  # every cell tied at 1: more than 2 k cells, so none is hot
  p <- cbind(x = 0:9 + 0.5, y = rep(0.5, 10))
  h <- hsper(p, prd = rep("a", 10), cell = 1, top = 0.1, org = c(0, 0))
  expect_identical(h$periods$hot, 0)
  expect_identical(h$periods$tiedout, 10)
  h2 <- hsper(p, prd = factor(rep(c("b", "a"), 5), levels = c("b", "a")),
              cell = 1, top = 0.1)
  expect_identical(h2$periods$period, c("b", "a"))
  expect_true(is.na(h2$periods$jacc[2]))
  set.seed(72)
  q <- cbind(x = runif(80, 0, 3000), y = runif(80, 0, 3000))
  ll <- eqinv(q, lon0 = -79.4, lat0 = 43.7)
  a <- hsper(ll, prd = rep(1:2, 40), cell = 500, lonlat = TRUE)
  b <- hsper(eqproj(ll), prd = rep(1:2, 40), cell = 500)
  expect_identical(a$hot, b$hot)
  # the centres are the planar ones, back about the points' own centre
  rad <- pi / 180
  c0 <- atan2(mean(sin(ll[, 1] * rad)), mean(cos(ll[, 1] * rad))) / rad
  back <- eqinv(as.matrix(b$cells[, c("x", "y")]), lon0 = c0,
                lat0 = mean(ll[, 2]))
  expect_equal(unname(as.matrix(a$cells[, c("x", "y")])),
               unname(as.matrix(back)),
               tolerance = 1e-12)
})

test_that("hsper() refuses bad input in words", {
  p <- cbind(x = 1:3, y = 1:3)
  errs("hsper", list(
    list("`cell`", p, prd = 1:3, cell = 0),
    list("`top` must be one share", p, prd = 1:3, cell = 1, top = 0),
    list("`prd` must have one label", p, prd = 1:2, cell = 1),
    list("`org`", p, prd = 1:3, cell = 1, org = 1),
    list("too many for the grid", p, prd = 1:3, cell = 1, org = c(-1e20, 0))
  ))
})
