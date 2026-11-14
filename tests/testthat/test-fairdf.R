# fairdf(), tmpaud() and lipaud(): against formulas written here from the
# papers (Wilson intervals from prop.test()), against brute force, against
# the recorded values of the owner's earlier package (dev/ref-fair.R), the
# disagreements with it, edge cases and refusals.

ref <- readRDS(test_path("ref-fair.rds"))
# the inputs of dev/ref-fair.R
fairin <- function() {
  set.seed(2)
  n <- 300
  g <- sample(c("A", "B", "C"), n, TRUE)
  y <- rbinom(n, 1, 0.4)
  p <- rbinom(n, 1, c(A = 0.5, B = 0.35, C = 0.2)[g])
  per <- sample(c("p1", "p2", "p3"), n, TRUE)
  cty <- sample(c("X", "Y"), n, TRUE)
  list(g = g, y = y, p = p, per = per, cty = cty)
}
d <- fairin()
z <- qnorm(0.975)
wil <- function(x, m) {
  if (m == 0) return(c(NaN, NaN))
  suppressWarnings(prop.test(x, m, correct = FALSE)$conf.int[1:2])
}
# Newcombe's method 10 from the two Wilson intervals
newc <- function(x1, m1, x2, m2) {
  a <- wil(x1, m1)
  b <- wil(x2, m2)
  p1 <- x1 / m1
  p2 <- x2 / m2
  c(p1 - p2 - sqrt((p1 - a[1])^2 + (b[2] - p2)^2),
    p1 - p2 + sqrt((a[2] - p1)^2 + (p2 - b[1])^2))
}
# the Gini coefficient as the mean absolute difference over twice the mean
gin <- function(x) sum(abs(outer(x, x, "-"))) / (2 * length(x)^2 * mean(x))

test_that("shares, gaps and ratios follow their formulas", {
  r <- fairdf(d$p, d$g, d$y, base = "A")
  t <- r$tab
  expect_identical(t$group, c("A", "B", "C"))
  expect_identical(r$base, "A")
  n <- as.vector(table(d$g))
  x <- as.vector(tapply(d$p, d$g, sum))
  expect_identical(t$n, as.double(n))
  expect_identical(t$npos, as.double(x))
  sh <- x / n
  expect_equal(t$ppos, sh, tolerance = 1e-15)
  expect_equal(t$dp, sh - sh[1], tolerance = 1e-15)
  expect_equal(t$di, sh / sh[1], tolerance = 1e-15)
  for (j in 1:3) {
    expect_equal(c(t$lo[j], t$hi[j]), wil(x[j], n[j]), tolerance = 1e-13)
  }
  for (j in 2:3) {
    expect_equal(c(t$dplo[j], t$dphi[j]), newc(x[j], n[j], x[1], n[1]),
                 tolerance = 1e-13)
    se <- sqrt(1 / x[j] - 1 / n[j] + 1 / x[1] - 1 / n[1])
    expect_equal(c(t$dilo[j], t$dihi[j]), sh[j] / sh[1] * exp(c(-z, z) * se),
                 tolerance = 1e-13)
  }
  expect_true(all(is.nan(unlist(t[1, c("dplo", "dphi", "dilo", "dihi")]))))
  # the rates given the outcomes
  tp <- as.vector(tapply(d$p[d$y == 1], d$g[d$y == 1], sum))
  np <- as.vector(table(d$g[d$y == 1]))
  fp <- as.vector(tapply(d$p[d$y == 0], d$g[d$y == 0], sum))
  nn <- as.vector(table(d$g[d$y == 0]))
  expect_equal(t$tpr, tp / np, tolerance = 1e-15)
  expect_equal(t$fpr, fp / nn, tolerance = 1e-15)
  expect_equal(t$dtpr, tp / np - tp[1] / np[1], tolerance = 1e-15)
  expect_equal(t$aod, (t$dtpr + t$dfpr) / 2, tolerance = 1e-15)
  for (j in 2:3) {
    a <- newc(tp[j], np[j], tp[1], np[1])
    b <- newc(fp[j], nn[j], fp[1], nn[1])
    expect_equal(c(t$dtprlo[j], t$dtprhi[j]), a, tolerance = 1e-13)
    expect_equal(c(t$dfprlo[j], t$dfprhi[j]), b, tolerance = 1e-13)
    # MOVER for the sum, halved
    s <- t$dtpr[j] + t$dfpr[j]
    mv <- c(s - sqrt((t$dtpr[j] - a[1])^2 + (t$dfpr[j] - b[1])^2),
            s + sqrt((a[2] - t$dtpr[j])^2 + (b[2] - t$dfpr[j])^2)) / 2
    expect_equal(c(t$aodlo[j], t$aodhi[j]), mv, tolerance = 1e-13)
    e <- c(t$dtpr[j], t$dfpr[j])
    expect_identical(t$eo[j], e[which.max(abs(e))])
  }
  s <- r$sum
  expect_named(s, c("di", "dp", "gini", "bas", "ginilo", "ginihi", "baslo",
                    "bashi", "aod", "eo"))
  expect_identical(s[["di"]], min(t$di[-1]))
  expect_identical(s[["dp"]], t$dp[-1][which.max(abs(t$dp[-1]))])
  expect_equal(s[["gini"]], gin(sh), tolerance = 1e-14)
  expect_equal(s[["bas"]], s[["dp"]] * s[["gini"]], tolerance = 1e-15)
  expect_true(all(is.nan(s[5:8])))
})

test_that("the values agree with the earlier package's", {
  r <- fairdf(d$p, d$g, d$y, base = "A")
  t <- r$tab
  expect_equal(t$di, unname(unlist(ref$di)[t$group]), tolerance = 1e-14)
  expect_equal(t$dp, unname(unlist(ref$dp)[t$group]), tolerance = 1e-14)
  expect_equal(r$sum[c("bas", "dp", "gini")], ref$bas[c("bas", "gap", "gini")],
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(t$dtpr, unname(ref$eo$tpr[t$group]), tolerance = 1e-14)
  expect_equal(t$dfpr, unname(ref$eo$fpr[t$group]), tolerance = 1e-14)
  expect_equal(r$sum[["eo"]], ref$eo$value, tolerance = 1e-14)
  expect_equal(t$aod, unname(ref$aod$aod[t$group]), tolerance = 1e-14)
  expect_equal(r$sum[["aod"]], ref$aod$value, tolerance = 1e-14)
  # the ratios of rmorie::morie_fairness_disparate_impact(d$p, d$g,
  # privileged = "A"), stored in ref/fairdf.txt by dev/ref/fairdf.R
  di <- refget("fairdf")
  expect_equal(c(sum(d$y), sum(d$p)), di$inp, tolerance = 1e-12,
               info = "inputs changed")
  expect_equal(t$di, unname(di$di[t$group]), tolerance = 1e-14)
})

test_that("disagreements with the earlier package are documented", {
  # 1. With no reference given, the earlier package compares equalized odds
  # and the average odds difference with the group with the highest true
  # positive rate, the ratio with the one with the highest share flagged:
  # one call, two references. fairdf() uses one base for every column.
  expect_identical(unname(ref$aodnull), "B")
  expect_identical(fairdf(d$p, d$g, d$y)$base, "A")
  # 2. A base group with no case flagged: the earlier package reports every
  # ratio as NaN; fairdf() divides as IEEE does, Inf over a positive share.
  expect_true(is.nan(ref$di0$b))
  r <- fairdf(c(0, 0, 1, 0), c("a", "a", "b", "b"), base = "a")
  expect_identical(r$tab$di, c(NaN, Inf))
  expect_true(is.nan(r$sum[["di"]]))
})

test_that("edge cases: one case, empty denominators, a constant decision", {
  r <- fairdf(c(1, 0, 1, 1), c("a", "b", "b", "b"), y = c(1, 0, 0, 1))
  # base by default: the largest share flagged, a (1 of 1)
  expect_identical(r$base, "a")
  t <- r$tab
  expect_identical(t$hi[1], 1)
  # no negatives in a: its false positive rate and every gap in it are NaN
  expect_true(is.nan(t$fpr[1]))
  expect_true(is.nan(t$dfpr[2]))
  expect_identical(t$eo[2], t$dtpr[2])
  expect_identical(t$eo[1], 0)
  # everyone flagged: shares 1, Gini 0, no amplification
  r <- fairdf(rep(TRUE, 6), rep(1:3, 2))
  expect_identical(r$sum[["gini"]], 0)
  expect_identical(r$sum[["bas"]], 0)
  expect_identical(r$tab$group, c("1", "2", "3"))
  expect_false("aod" %in% names(r$sum))
  expect_length(r$tab, 12L)
})

test_that("the bootstrap is reproducible and brackets the estimates", {
  set.seed(5)
  g <- rep(c("a", "b", "c"), c(150, 120, 90))
  p <- rbinom(360, 1, c(a = 0.5, b = 0.3, c = 0.2)[g])
  a <- fairdf(p, g, nboot = 499, seed = 3)$sum
  expect_identical(fairdf(p, g, nboot = 499, seed = 3)$sum, a)
  expect_false(identical(fairdf(p, g, nboot = 499, seed = 4)$sum, a))
  set.seed(9)
  b1 <- fairdf(p, g, nboot = 99)$sum
  set.seed(9)
  expect_identical(fairdf(p, g, nboot = 99)$sum, b1)
  expect_lt(a[["ginilo"]], a[["gini"]])
  expect_gt(a[["ginihi"]], a[["gini"]])
  expect_lt(a[["baslo"]], a[["bas"]])
  expect_gt(a[["bashi"]], a[["bas"]])
  # coverage of the Gini interval of the generating shares, over replicate
  # samples: a percentile interval near its level (here 95%) for n this
  # large; 40 samples, at least 30 covering is a loose check
  tru <- gin(c(0.5, 0.3, 0.2))
  cov <- vapply(1:40, function(i) {
    p <- rbinom(360, 1, c(a = 0.5, b = 0.3, c = 0.2)[g])
    s <- fairdf(p, g, nboot = 199, seed = i)$sum
    s[["ginilo"]] <= tru && tru <= s[["ginihi"]]
  }, logical(1))
  expect_gt(sum(cov), 29)
})

test_that("tmpaud() reproduces the earlier package's temporal audit", {
  r <- tmpaud(d$p, d$g, d$per, d$cty, base = "A")
  cl <- r$cells
  rc <- ref$tmp$cells
  rc <- rc[order(rc$city, rc$period), ]
  expect_identical(cl$city, rc$city)
  expect_identical(cl$period, rc$period)
  expect_identical(cl$n, rc$n)
  expect_equal(cl$di, rc$dir, tolerance = 1e-14)
  expect_equal(cl$dp, rc$pg, tolerance = 1e-14)
  expect_equal(cl$gini, rc$gini, tolerance = 1e-14)
  expect_equal(cl$bas, rc$bas, tolerance = 1e-14)
  cy <- r$cities
  pc <- ref$tmp$city[cy$city]
  expect_equal(cy$di, vapply(pc, `[[`, 0, "mean_dir"), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(cy$dirng, vapply(pc, `[[`, 0, "dir_range"), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_identical(cy$ngt1, vapply(pc, `[[`, 0L, "periods_dir_gt1"),
                   ignore_attr = TRUE)
  expect_equal(r$sum[["rng"]], ref$tmp$rng, tolerance = 1e-14)
  expect_equal(r$sum[["spread"]], ref$tmp$spread, tolerance = 1e-14)
  # each cell is fairdf() on its cases
  k <- d$cty == "Y" & d$per == "p2"
  f <- fairdf(d$p[k], d$g[k], base = "A")$sum
  expect_equal(unlist(cl[cl$city == "Y" & cl$period == "p2",
                         c("di", "dp", "gini", "bas")]),
               f[c("di", "dp", "gini", "bas")], tolerance = 1e-14,
               ignore_attr = TRUE)
  # one city: the earlier package reports a spread of 0 between cities;
  # there is no spread to report, and tmpaud() says so with NaN
  expect_identical(ref$spread1, 0)
  r1 <- tmpaud(d$p, d$g, d$per, base = "A")
  expect_identical(r1$cities$city, "all")
  expect_true(is.nan(r1$sum[["spread"]]))
})

test_that("tmpaud() skips the cells it cannot audit", {
  p <- c(1, 0, 1, 1, 0, 1, 0, 0)
  g <- c("a", "b", "a", "a", "b", "b", "a", "b")
  per <- c(1, 1, 2, 2, 3, 3, 3, 3)
  r <- tmpaud(p, g, per)
  # the base by the pooled shares: a, 3 of 4
  expect_identical(r$base, "a")
  # period 2 holds group a only
  expect_identical(r$cells$ok, c(TRUE, FALSE, TRUE))
  expect_true(all(is.nan(unlist(r$cells[2, c("di", "dp", "gini", "bas")]))))
  expect_identical(r$cities$ncell, 2L)
  # base b: period 2 lacks it, period 1 has b's share 0 (no finite ratio)
  r <- tmpaud(p, g, per, base = "b")
  expect_identical(r$cells$ok, c(TRUE, FALSE, TRUE))
  expect_true(is.nan(r$cells$di[1]))
  expect_identical(r$cells$di[3], 0)
  expect_identical(r$cities$ngt1, 0L)
})

test_that("lipaud() checks every pair as brute force does", {
  set.seed(4)
  x <- matrix(rnorm(30), 10)
  h <- plogis(x %*% c(1, -1, 0.5))[, 1]
  dd <- dist(x)
  r <- lipaud(h, dd, lip = 0.3)
  m <- as.matrix(dd)
  ij <- which(lower.tri(m), arr.ind = TRUE)
  ij <- ij[order(ij[, 2], ij[, 1]), ]
  gap <- abs(h[ij[, 1]] - h[ij[, 2]])
  dv <- m[ij]
  expect_identical(r$pairs$i, as.integer(ij[, 2]))
  expect_identical(r$pairs$j, as.integer(ij[, 1]))
  expect_equal(r$pairs$gap, gap, tolerance = 1e-15)
  expect_equal(r$lreq, max(gap / dv), tolerance = 1e-15)
  expect_identical(r$nviol, sum(gap > 0.3 * dv))
  expect_equal(r$rate, mean(gap > 0.3 * dv), tolerance = 1e-15)
  expect_identical(r$ok, r$nviol == 0L)
  expect_identical(r$maxgap, max(gap))
  w <- which.max(gap / dv)
  expect_identical(r$worst, c(r$pairs$i[w], r$pairs$j[w]))
  # a symmetric matrix and an explicit list of pairs give the same audit
  expect_identical(lipaud(h, m, lip = 0.3), r)
  pr <- cbind(r$pairs$i, r$pairs$j)
  expect_identical(lipaud(h, r$pairs$d, pr, lip = 0.3), r)
  # at the constant it reports, nothing breaks the rule
  expect_true(lipaud(h, dd, lip = r$lreq * (1 + 1e-12))$ok)
})

test_that("lipaud(): distance 0, and the earlier package's indexing", {
  h <- c(1, 2, 3, 4, 5, 6, 7, 8)
  # the earlier package numbers people from 0: its pairs 0-4, 1-5 and 2-7
  # are 1-5, 2-6 and 3-8 here; at distance 0, a gap breaks the rule for
  # every constant
  pr <- matrix(c(1, 2, 3, 5, 6, 8), ncol = 2)
  r <- lipaud(h, c(1, 2, 0), pr, lip = 1.5)
  expect_identical(r$lreq, ref$lip$L_required)
  expect_identical(r$nviol, ref$lip$n_violations)
  expect_identical(r$maxgap, ref$lip$max_gap)
  expect_identical(r$pairs$i[which.max(r$pairs$gap)] - 1L, ref$lip$max_pair_i)
  expect_identical(r$pairs$ratio, c(4, 2, Inf))
  expect_identical(r$worst, c(3L, 8L))
  # two people alike and scored alike comply
  pr <- matrix(c(1, 1, 2, 3), 2)
  r <- lipaud(c(0.2, 0.2, 0.5), c(0, 1), pr, lip = 0)
  expect_identical(r$pairs$ratio, c(0, 0.3))
  expect_identical(r$pairs$viol, c(FALSE, TRUE))
  expect_identical(lipaud(1:3 / 4, pairs = matrix(1:2, 1))$pairs$d, 1)
})

test_that("fairdf(), tmpaud() and lipaud() refuse what they cannot use", {
  errs("fairdf", list(
    list("at least one case", numeric(0), character(0)),
    list("must be 0 or 1", c(0, 2), c("a", "b")),
    list("at least two groups", c(0, 1), c("a", "a")),
    list("one label per value", c(0, 1), c("a", "b", "c")),
    list("must be 0 or 1", c(0, 1), c("a", "b"), y = c(1, 3)),
    list("`base` names no group", c(0, 1), c("a", "b"), base = "c"),
    list("strictly between", c(0, 1), c("a", "b"), conf = 2),
    list("whole number", c(0, 1), c("a", "b"), nboot = -1),
    list("Decode it", structure(c(0, 1), labels = c(no = 0, yes = 1)),
         c("a", "b")),
    list("not a matrix", matrix(c(0, 1, 1, 0), 2), 1:4)
  ))
  errs("tmpaud", list(
    list("at least one case", numeric(0), character(0), character(0)),
    list("at least two groups", c(0, 1), c("a", "a"), 1:2),
    list("one label per value", c(0, 1), c("a", "b"), 1),
    list("one label per value", c(0, 1), c("a", "b"), 1:2, city = 1)
  ))
  errs("lipaud", list(
    list("give `d`", 1:3),
    list("square matrix", 1:3, diag(2)),
    list("square matrix", 1:3, dist(1:4)),
    list("symmetric", 1:2, matrix(c(0, 1, 2, 0), 2)),
    list("finite distances", 1:2, matrix(c(0, -1, -1, 0), 2)),
    list("at least two cases", 1, matrix(0, 1, 1)),
    list("two columns", 1:3, NULL, matrix(1:3, 1)),
    list("case numbers from 1 to 3", 1:3, NULL, matrix(c(1, 4), 1)),
    list("must have length 1", 1:3, 1:2, matrix(c(1, 2), 1)),
    list("`lip` must be at least 0", 1:3, dist(1:3), lip = -1),
    list("finite numbers", c(1, NA), dist(1:2))
  ))
})
