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
  expect_gt(r$thr$se[1], 0)
  expect_lt(r$thr$lo[1], r$thr$thr[1])
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

# Counts made exactly from the model of Anwar and Fang (2006, Section 2):
# among people of group r, a share pi carry contraband; an officer of
# group o sees a signal N(1, 1) from the guilty and N(0, 1) from the rest
# and searches when the chance of guilt passes t[r, o].
afsim <- function(t, n = 1e5) {
  pi <- c(0.2, 0.3)
  d <- expand.grid(off = seq_len(ncol(t)), grp = seq_len(nrow(t)))
  ph <- pi[d$grp]
  tt <- t[cbind(d$grp, d$off)]
  x <- (qlogis(tt) - qlogis(ph)) + 0.5
  A <- ph * pnorm(x - 1, lower.tail = FALSE)
  B <- (1 - ph) * pnorm(x, lower.tail = FALSE)
  data.frame(d, stops = n, srch = n * (A + B), hits = n * A)
}

test_that("the joint rank-order test holds without prejudice only", {
  # thresholds that depend on the officer group only: no prejudice
  t0 <- rbind(c(0.2, 0.3), c(0.2, 0.3))
  d <- afsim(t0)
  r <- thrtst(d$stops, d$srch, d$hits, d$grp, off = d$off, rank = "joint")
  expect_false(any(r$rank$rej))
  # officers of group 1 search more often and find less, in both groups
  expect_identical(r$rank$rank, c(1, 2, 2, 1))
  expect_named(r$rank, c("meas", "off", "rank", "rej"))
  # officer group 2 is the stricter with people of group 1 and the more
  # lenient with people of group 2: prejudice, the rankings flip
  t1 <- rbind(c(0.2, 0.3), c(0.3, 0.2))
  d1 <- afsim(t1)
  r1 <- thrtst(d1$stops, d1$srch, d1$hits, d1$grp, off = d1$off, rank = "joint")
  expect_true(all(r1$rank$rej))
  expect_true(all(is.na(r1$rank$rank)))
  # the pairwise form agrees here
  p1 <- thrtst(d1$stops, d1$srch, d1$hits, d1$grp, off = d1$off)
  expect_true(all(p1$rank$rej))
  # the test of monolithic behaviour: Pearson's chi-square by group
  m <- r1$mono
  k <- m$grp == "1" & m$meas == "search"
  ct <- suppressWarnings(chisq.test(
    cbind(d1$srch[d1$grp == 1], (d1$stops - d1$srch)[d1$grp == 1]),
    correct = FALSE
  ))
  expect_equal(m$chisq[k], unname(ct$statistic), tolerance = 1e-12)
  expect_identical(m$df[k], 1)
  expect_null(thrtst(d$stops, d$srch, d$hits, d$grp)$mono)
})

test_that("the joint test finds a cycle that no pair shows", {
  # search rates of officer groups a, b, c among people of groups 1, 2, 3:
  # a > b in group 1, b > c in group 2, c > a in group 3, the rest within
  # chance
  sr <- c(72, 48, 60, 60, 72, 48, 48, 60, 72)
  g <- rep(1:3, each = 3)
  o <- rep(c("a", "b", "c"), 3)
  r <- thrtst(rep(600, 9), sr, sr / 2, g, off = o, rank = "joint")
  expect_true(all(r$rank$rej[r$rank$meas == "search"]))
  expect_false(any(r$rank$rej[r$rank$meas == "hit"]))
  # the one-sided z of equation 7 for a against b in group 1
  pa <- 72 / 600
  pb <- 48 / 600
  z <- (pa - pb) / sqrt(pa * (1 - pa) / 599 + pb * (1 - pb) / 599)
  expect_gt(z, qnorm(0.95))
  p <- thrtst(rep(600, 9), sr, sr / 2, g, off = o)
  expect_false(any(p$rank$rej))
})

# Three groups with the same counts in five areas: the posterior is the
# same for every group, so the differences of their thresholds are 0 but
# for the Monte Carlo error.
bsym <- function() {
  d <- expand.grid(g = c("a", "b", "c"), a = 1:5, stringsAsFactors = FALSE)
  d$stops <- c(60, 50, 40, 55, 45)[d$a]
  d$srch <- c(12, 9, 8, 11, 7)[d$a]
  d$hits <- c(5, 3, 3, 4, 2)[d$a]
  d
}

test_that("the Bayesian fit is symmetric, reproducible and reports its draws", {
  d <- bsym()
  b <- thrtst(d$stops, d$srch, d$hits, d$g, area = d$a, fit = "bayes",
              iter = 500, chains = 2, seed = 3)
  expect_named(
    b$thr,
    c("group", "thr", "se", "lo", "hi", "diff", "dse", "dlo", "dhi")
  )
  expect_null(b$lrt)
  expect_identical(b$fit$model, "bayes")
  expect_identical(b$fit$chains, 2L)
  k <- b$post$param %in% c("diff[b]", "diff[c]")
  mc <- b$post$sd[k] / sqrt(b$post$ess[k])
  expect_true(all(abs(b$post$mean[k]) < 4 * mc))
  expect_lt(max(b$post$rhat, na.rm = TRUE), 1.05)
  expect_gte(min(b$post$ess[1:6], na.rm = TRUE), 100)
  expect_true(all(b$thr$lo < b$thr$thr))
  expect_true(all(b$thr$thr < b$thr$hi))
  expect_equal(b$thr$thr, b$post$mean[1:3], tolerance = 1e-15)
  # 3 groups, 4 hyperparameters, 4 areas but the pinned one, 15 cells
  expect_identical(nrow(b$post), 6L * 3L + 4L + 2L * 4L + 15L)
  expect_identical(nrow(b$cells), 15L)
  expect_identical(
    b,
    thrtst(d$stops, d$srch, d$hits, d$g, area = d$a, fit = "bayes",
           iter = 500, chains = 2, seed = 3)
  )
})

test_that("the Bayesian fit sees a lower threshold", {
  d <- bsym()
  # group c searched twice as often, with as many hits
  sr <- ifelse(d$g == "c", 2 * d$srch, d$srch)
  b <- thrtst(d$stops, sr, d$hits, d$g, area = d$a, fit = "bayes",
              iter = 500, seed = 3)
  dc <- b$thr$diff[3]
  expect_lt(dc, 0)
  expect_lt(dc, b$thr$diff[2])
})

test_that("the Bayesian fit refuses unidentified data and unmixed chains", {
  d <- bsym()
  errs(
    "thrtst",
    list(
      list(
        "needs more cells with stops than 2",
        d$stops[d$g != "c"], d$srch[d$g != "c"], d$hits[d$g != "c"],
        d$g[d$g != "c"], area = d$a[d$g != "c"], fit = "bayes", seed = 1
      ),
      list(
        "the chains have not mixed",
        d$stops, d$srch, d$hits, d$g, area = d$a, fit = "bayes",
        iter = 20, seed = 1
      ),
      list("`fit` must be one of", c(10, 10), c(2, 3), c(1, 1), c(1, 2),
           fit = "mcmc"),
      list("`chains`", c(10, 10), c(2, 3), c(1, 1), c(1, 2), chains = 1),
      list("`iter`", c(10, 10), c(2, 3), c(1, 1), c(1, 2), iter = 10)
    )
  )
})

test_that("on data with one threshold per group the Bayesian fit sits below", {
  # the behaviour the help page warns of: counts made with one threshold
  # per group and area offsets (three groups, six areas); the model with a
  # threshold for every cell puts every group's threshold below the
  # maximum-likelihood one
  d <- expand.grid(g = 1:3, a = 1:6)
  ag <- c(-1.2, -0.6, -0.9)
  bg <- c(0.4, 0.2, 0.3)
  cg <- c(-1, -1.4, -0.8)
  aa <- c(0, 0.3, -0.4, 0.2, -0.1, 0.25)
  ba <- c(0, 0.15, -0.1, 0.05, 0.1, -0.05)
  ca <- c(0, 0.2, -0.3, 0.1, -0.15, 0.05)
  ph <- plogis(ag[d$g] + aa[d$a])
  dl <- exp(bg[d$g] + ba[d$a])
  gg <- cg[d$g] + ca[d$a]
  x <- (gg - qlogis(ph)) / dl + dl / 2
  A <- ph * pnorm(x - dl, lower.tail = FALSE)
  B <- (1 - ph) * pnorm(x, lower.tail = FALSE)
  st <- ifelse(d$a == 1, 1000, 500)
  sr <- round(st * (A + B))
  ht <- round(st * A)
  m <- thrtst(st, sr, ht, d$g, area = d$a)
  b <- thrtst(st, sr, ht, d$g, area = d$a, fit = "bayes", iter = 500,
              chains = 2, seed = 1)
  expect_true(all(b$thr$thr < m$thr$thr))
})
