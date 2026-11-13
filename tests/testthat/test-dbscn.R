# DBSCAN. The kernel is checked by hand in src/t_hot.cpp; here the binding
# against dbscan::dbscan(), stored.

test_that("clusters are those of dbscan::dbscan()", {
  set.seed(62)
  q <- rbind(cbind(rnorm(15, 0, 0.3), rnorm(15, 0, 0.3)),
             cbind(rnorm(15, 2, 0.3), rnorm(15, 1, 0.3)),
             cbind(runif(10, -1, 3), runif(10, -1, 2)))
  # dbscan 1.2.5: dbscan::dbscan(q, eps = 0.4, minPts = 4)$cluster and
  # with borderPoints = FALSE
  ref <- c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
           2L, 2L, 2L, 0L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L,
           0L, 1L, 0L, 0L, 0L, 0L, 2L, 1L, 1L)
  nob <- c(1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L,
           0L, 2L, 2L, 0L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 0L,
           0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L)
  d <- dbscn(q, eps = 0.4, minpts = 4)
  expect_s3_class(d, "dbscn")
  expect_named(d, c("cl", "core", "sizes"))
  expect_identical(d$cl, ifelse(ref == 0L, NA_integer_, ref))
  expect_identical(d$sizes, tabulate(ref[ref > 0]))
  b <- dbscn(q, eps = 0.4, minpts = 4, border = FALSE)
  expect_identical(b$cl, ifelse(nob == 0L, NA_integer_, nob))
  # core points: at least minpts within eps, themselves included
  dm <- as.matrix(dist(q))
  expect_identical(d$core, unname(rowSums(dm <= 0.4) >= 4))
  expect_output(print(d), "dbscn: 2 clusters of 40 points, 7 noise")
  expect_output(print(dbscn(q, eps = 0.01, minpts = 2)), "0 clusters")
  ll <- cbind(lon = -79 + q[, 1] / 100, lat = 43 + q[, 2] / 100)
  expect_identical(dbscn(ll, eps = 300, minpts = 4, lonlat = TRUE),
                   dbscn(eqproj(ll), eps = 300, minpts = 4))
})

test_that("dbscn() refuses bad input in words", {
  p <- cbind(x = 1:3, y = 1:3)
  errs("dbscn", list(
    list("`eps`", p, eps = 0),
    list("`minpts`", p, eps = 1, minpts = 0),
    list("`border`", p, eps = 1, border = NA)
  ))
})
