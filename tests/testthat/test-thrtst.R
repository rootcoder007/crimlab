# Tests of discrimination in searches. The kernel is checked by hand in
# src/t_thr.cpp; here the binding against base R's arithmetic and counts
# made exactly from the threshold model.

test_that("hit rates and the chi-square follow the formulas", {
  st <- c(100, 100, 80)
  sr <- c(20, 40, 16)
  ht <- c(8, 10, 4)
  r <- thrtst(st, sr, ht, c("a", "b", "c"))
  expect_s3_class(r, "thrtst")
  expect_equal(r$hit$hrate, ht / sr, tolerance = 1e-15)
  expect_equal(r$hit$srate, sr / st, tolerance = 1e-15)
  h <- ht / sr
  expect_equal(r$hit$diff, h - h[1], tolerance = 1e-14)
  se <- sqrt(h * (1 - h) / sr + h[1] * (1 - h[1]) / sr[1])
  expect_equal(r$hit$se[-1], se[-1], tolerance = 1e-14)
  expect_equal(
    r$hit$p[-1],
    2 * pnorm(-abs((h - h[1]) / se))[-1],
    tolerance = 1e-12
  )
  expect_true(is.na(r$hit$z[1]))
  ct <- suppressWarnings(chisq.test(cbind(ht, sr - ht), correct = FALSE))
  expect_equal(
    unname(r$kpt[["chisq"]]),
    unname(ct$statistic),
    tolerance = 1e-13
  )
  expect_equal(unname(r$kpt[["p"]]), ct$p.value, tolerance = 1e-12)
  expect_null(r$thr)
  expect_output(print(r), "not identified")
  r2 <- thrtst(st, sr, ht, c("a", "b", "c"), ref = "b")
  expect_identical(r2$ref, "b")
})

# counts that are exactly the model's expected search and hit counts
thrsim <- function(dt) {
  ag <- c(-1.2, -0.6)
  bg <- c(0.4, 0.2)
  cg <- c(-1, -1 + dt)
  aa <- c(0, 0.3, -0.4, 0.2)
  ba <- c(0, 0.15, -0.1, 0.05)
  ca <- c(0, 0.2, -0.3, 0.1)
  d <- expand.grid(g = 1:2, a = 1:4)
  ph <- plogis(ag[d$g] + aa[d$a])
  dl <- exp(bg[d$g] + ba[d$a])
  gg <- cg[d$g] + ca[d$a]
  x <- (gg - qlogis(ph)) / dl + dl / 2
  A <- ph * pnorm(x - dl, lower.tail = FALSE)
  B <- (1 - ph) * pnorm(x, lower.tail = FALSE)
  list(
    d = d,
    stops = rep(20000, 8),
    srch = 20000 * (A + B),
    hits = 20000 * A,
    thr = tapply(plogis(gg), d$g, mean)
  )
}

test_that("the threshold model gives back the thresholds of the counts", {
  s <- thrsim(0.5)
  r <- thrtst(s$stops, s$srch, s$hits, s$d$g, area = s$d$a)
  expect_equal(r$thr$thr, unname(as.vector(s$thr)), tolerance = 1e-6)
  expect_equal(r$cells$srate, s$srch / s$stops, tolerance = 1e-6)
  expect_true(all(r$fit$conv))
  expect_gt(r$lrt[["stat"]], 0)
  expect_identical(r$lrt[["cells"]], 8)
  expect_true(r$thr$se[1] > 0 && r$thr$lo[1] < r$thr$thr[1])
  s0 <- thrsim(0)
  r0 <- thrtst(s0$stops, s0$srch, s0$hits, s0$d$g, area = s0$d$a)
  expect_lt(r0$lrt[["stat"]], 1e-4)
  expect_output(print(r0), "likelihood ratio")
})

test_that("the rank-order test reads officer groups whose ranking flips", {
  r <- thrtst(
    c(100, 100, 100, 100),
    c(10, 20, 30, 10),
    c(5, 5, 5, 5),
    grp = c(1, 1, 2, 2),
    off = c("x", "y", "x", "y")
  )
  k <- r$rank
  expect_identical(k$meas, c("search", "hit"))
  expect_equal(k$zmin[1], -2, tolerance = 1e-14)
  expect_equal(k$zmax[1], 0.2 / sqrt(0.003), tolerance = 1e-14)
  expect_identical(k$rej, c(TRUE, FALSE))
  expect_output(print(r), "rank order")
})

test_that("thrtst() refuses bad input in words", {
  errs(
    "thrtst",
    list(
      list("`srch` must be at most", c(10, 10), c(11, 1), c(1, 1), c(1, 2)),
      list(
        "`grp` must hold at least two",
        c(10, 10),
        c(1, 1),
        c(1, 1),
        c(1, 1)
      ),
      list("every group needs searches", c(10, 10), c(0, 1), c(0, 1), c(1, 2)),
      list("`level`", c(10, 10), c(1, 1), c(1, 1), c(1, 2), level = 1),
      list(
        "`ref` names no group",
        c(10, 10),
        c(1, 1),
        c(1, 1),
        c(1, 2),
        ref = "z"
      ),
      list(
        "`off` must hold",
        c(10, 10),
        c(1, 1),
        c(1, 1),
        c(1, 2),
        off = c(1, 1)
      )
    )
  )
})
