# Every arrangement of the treated within blocks, enumerated with combn():
# the values of the difference in means (or the block-weighted
# difference) over them.
bp_enum <- function(y, t, b, strat = FALSE) {
  idx <- split(seq_along(y), factor(b, unique(b)))
  cb <- lapply(idx, function(i) {
    k <- sum(t[i])
    if (k == 0L) return(matrix(integer(0), 0L, 1L))
    if (length(i) == 1L || k == length(i)) return(matrix(i, k, 1L))
    utils::combn(i, k)
  })
  g <- expand.grid(lapply(cb, function(m) seq_len(ncol(m))))
  apply(g, 1L, function(r) {
    tt <- integer(length(y))
    for (j in seq_along(cb)) tt[cb[[j]][, r[j]]] <- 1L
    if (!strat) return(mean(y[tt == 1L]) - mean(y[tt == 0L]))
    sum(vapply(idx, function(i) {
      u <- tt[i]
      length(i) / length(y) * (mean(y[i][u == 1L]) - mean(y[i][u == 0L]))
    }, 0))
  })
}

bp_dat <- data.frame(
  y = c(3.1, 0.4, 2.2, 5.0, 1.7, 2.9, 0.8, 4.4, 1.1, 3.6),
  t = c(1, 0, 0, 1, 1, 0, 1, 0, 0, 1),
  b = rep(c("a", "b", "c"), c(4, 4, 2))
)

test_that("blkprm() with exact gives the enumerated permutation p-value", {
  d <- bp_dat
  obs <- mean(d$y[d$t == 1]) - mean(d$y[d$t == 0])
  st <- bp_enum(d$y, d$t, d$b)
  tol <- 1e-9 * (abs(obs) + max(abs(d$y)))
  r <- blkprm(d, "y", "t", "b", exact = TRUE)
  expect_equal(r$est, obs, tolerance = 1e-14)
  expect_identical(r$nperm, 72)
  expect_equal(r$p, mean(abs(st) >= abs(obs) - tol),
               tolerance = 1e-15)
  expect_identical(r$nge, sum(abs(st) >= abs(obs) - tol) + 0)
  expect_equal(sort(r$dist), sort(st), tolerance = 1e-14)
  expect_identical(r$seed, NA_real_)
  g <- blkprm(d, "y", "t", "b", alt = "greater", exact = TRUE)
  expect_equal(g$p, mean(st >= obs - tol),
               tolerance = 1e-15)
  l <- blkprm(d, "y", "t", "b", alt = "less", exact = TRUE)
  expect_equal(l$p, mean(st <= obs + tol),
               tolerance = 1e-15)
})

test_that("blkprm() weights the blocks' differences for stat = 'strat'", {
  d <- bp_dat
  ib <- split(seq_len(10), d$b)
  obs <- sum(vapply(ib, function(i) {
    length(i) / 10 * (mean(d$y[i][d$t[i] == 1]) - mean(d$y[i][d$t[i] == 0]))
  }, 0))
  st <- bp_enum(d$y, d$t, d$b, strat = TRUE)
  tol <- 1e-9 * (abs(obs) + max(abs(d$y)))
  r <- blkprm(d, "y", "t", "b", stat = "strat", exact = TRUE)
  expect_equal(r$est, obs, tolerance = 1e-14)
  expect_equal(r$p, mean(abs(st) >= abs(obs) - tol),
               tolerance = 1e-15)
  expect_equal(sort(r$dist), sort(st), tolerance = 1e-14)
})

test_that("blkprm()'s drawn arrangements agree with the enumeration", {
  d <- bp_dat
  ex <- blkprm(d, "y", "t", "b", exact = TRUE)$p
  r <- blkprm(d, "y", "t", "b", nperm = 20000, seed = 9)
  expect_identical(r$nperm, 20000)
  expect_identical(r$p, (1 + r$nge) / 20001)
  expect_lt(abs(r$p - ex), 4 * sqrt(ex * (1 - ex) / 20000) + 1e-3)
  # a drawn arrangement keeps every block's treated count, so the
  # permuted values are among the enumerated ones
  st <- bp_enum(d$y, d$t, d$b)
  near <- vapply(r$dist[1:200], function(v) min(abs(st - v)), 0)
  expect_lt(max(near), 1e-12)
  expect_identical(blkprm(d, "y", "t", "b", nperm = 100, seed = 9)$dist,
                   r$dist[1:100])
  set.seed(3)
  a <- blkprm(d, "y", "t", "b", nperm = 50)
  expect_identical(blkprm(d, "y", "t", "b", nperm = 50, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_match(out, "^blkprm: difference in means ")
  expect_match(out, "\\(two.sided, 20000 arrangements within 3 blocks\\)$")
  expect_identical(v, r)
})

test_that("blkprm() refuses what the test cannot use", {
  d <- bp_dat
  expect_error(blkprm(d, "y", "y", "b", exact = TRUE),
               "blkprm: column `y` must hold 0 or 1.", fixed = TRUE)
  d1 <- transform(d, t = 1)
  expect_error(blkprm(d1, "y", "t", "b", exact = TRUE),
               "blkprm: column `t` must hold treated and control units.",
               fixed = TRUE)
  d2 <- d
  d2$t[9:10] <- c(1, 1)
  expect_error(blkprm(d2, "y", "t", "b", stat = "strat", exact = TRUE),
               paste("blkprm: with `stat = \"strat\"` every block needs",
                     "treated and control units."), fixed = TRUE)
  big <- data.frame(y = seq_len(40), t = rep(0:1, 20), b = 1)
  expect_error(blkprm(big, "y", "t", "b", exact = TRUE),
               paste("blkprm: `exact` enumerates at most 1000000",
                     "arrangements; these blocks have 137846528820."),
               fixed = TRUE)
  expect_error(blkprm(d, "y", "t", "b", nperm = 0, seed = 1),
               "blkprm: `nperm` must be between 1 and 1e+07.", fixed = TRUE)
  expect_error(blkprm(d, "y", "t", "b", stat = "mean"),
               paste("blkprm: `stat` must be one of \"diff\", \"strat\",",
                     "not \"mean\"."), fixed = TRUE)
  expect_error(blkprm(d, "y", "t", "b", exact = NA),
               "blkprm: `exact` must be TRUE or FALSE.", fixed = TRUE)
})
