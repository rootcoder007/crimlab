# The space-time permutation scan. The kernel is checked against brute
# force over every cylinder in src/t_kscn.cpp; here the binding against the
# scan written out in R.

scanr <- function(p, t, r, tw, minc = 2) {
  n <- nrow(p)
  d <- as.matrix(dist(p))
  best <- c(llr = 0)
  for (c in seq_len(n)) for (rr in r) for (s in unique(t)) for (l in tw) {
    ins <- d[c, ] <= rr
    inw <- t >= s & t < s + l
    k <- sum(ins & inw)
    e <- sum(ins) * sum(inw) / n
    if (k >= minc && k > e) {
      v <- k * log(k / e) + if (n > k) (n - k) * log((n - k) / (n - e)) else 0
      if (v > best[["llr"]]) best <- c(llr = v, n = k, exp = e, r = rr,
                                     start = s, end = s + l)
    }
  }
  best
}

test_that("the most likely cluster is the scan's largest ratio", {
  set.seed(81)
  n <- 60
  p <- cbind(x = runif(n, 0, 10), y = runif(n, 0, 10))
  t <- sample(0:29, n, replace = TRUE)
  s <- kscan(p, t = t, r = c(1, 2.5), tw = c(3, 8), nsim = 0)
  b <- scanr(p, t, c(1, 2.5), c(3, 8))
  expect_s3_class(s, "kscan")
  expect_named(s, c("cl", "memb", "n", "nsim"))
  expect_equal(s$cl$llr[1], b[["llr"]], tolerance = 1e-13)
  expect_identical(c(s$cl$n[1], s$cl$r[1], s$cl$start[1], s$cl$end[1]),
                   unname(b[c("n", "r", "start", "end")]))
  expect_equal(s$cl$exp[1], b[["exp"]], tolerance = 1e-15)
  expect_equal(s$cl$rr[1], (b[["n"]] / b[["exp"]]) /
                 ((n - b[["n"]]) / (n - b[["exp"]])), tolerance = 1e-15)
  expect_identical(sum(s$memb == 1, na.rm = TRUE), as.integer(b[["n"]]))
  expect_true(is.nan(s$cl$p[1]))
  expect_output(print(s), "kscan: 1 clusters of 60 incidents")
})

test_that("a planted outbreak is found and tested", {
  set.seed(82)
  p <- cbind(x = c(runif(150, 0, 20), rnorm(20, 5, 0.3)),
             y = c(runif(150, 0, 20), rnorm(20, 15, 0.3)))
  t <- as.Date("2024-01-01") +
    c(sample(0:99, 150, TRUE), sample(40:44, 20, TRUE))
  s <- kscan(p, t = t, r = c(1, 3), tw = c(5, 20), top = 3, nsim = 99,
             seed = 1)
  expect_identical(s, kscan(p, t = as.numeric(t), r = c(1, 3),
                            tw = c(5, 20), top = 3, nsim = 99, seed = 1))
  expect_lt(abs(s$cl$x[1] - 5), 1)
  expect_identical(s$cl$p[1], 0.01)
  expect_true(all(diff(s$cl$llr) <= 0))
  expect_true(all(s$memb[151:170] == 1, na.rm = TRUE))
  tm <- as.POSIXct(t) + 60
  expect_identical(kscan(p, t = tm, r = c(1, 3), tw = c(5, 20), nsim = 0)$cl$n,
                   kscan(p, t = t, r = c(1, 3), tw = c(5, 20), nsim = 0)$cl$n)
  ll <- eqinv(p * 100, lon0 = -79.4, lat0 = 43.7)
  a <- kscan(ll, t = t, r = c(100, 300), tw = c(5, 20), nsim = 0,
             lonlat = TRUE)
  expect_identical(a$memb, kscan(eqproj(ll), t = t, r = c(100, 300),
                                 tw = c(5, 20), nsim = 0)$memb)
})

test_that("kscan() refuses bad input in words", {
  p <- cbind(x = 1:4, y = c(2, 1, 4, 3))
  errs("kscan", list(
    list("at least 3 incidents", p[1:2, ], t = 1:2, r = 1, tw = 1),
    list("`r` must increase", p, t = 1:4, r = c(2, 1), tw = 1),
    list("`tw`", p, t = 1:4, r = 1, tw = 0),
    list("`top`", p, t = 1:4, r = 1, tw = 1, top = 0),
    list("`minc`", p, t = 1:4, r = 1, tw = 1, minc = 0)
  ))
})
