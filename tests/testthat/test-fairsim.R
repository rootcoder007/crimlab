# fairsim() and nordet(): the simulated records against the rules that
# make them, the audits recovering the injected disparity, the noisy-OR
# probabilities against brute force and the earlier package's recorded
# values, and refusals.

ref <- readRDS(test_path("ref-fair.rds"))

test_that("the records follow the rules that make them", {
  d <- fairsim(30000, groups = c("r", "s", "t"), props = c(2, 1, 1),
               nareas = 11, base = 0.4, bias = 0.3, seed = 5)
  expect_named(d, c("area", "group", "truth", "det", "risk"))
  expect_identical(nrow(d), 30000L)
  # segregated areas: area a holds group (a - 1) mod 3 + 1 only
  gi <- match(d$group, c("r", "s", "t"))
  expect_identical((d$area - 1L) %% 3L + 1L, gi)
  expect_true(all(d$area >= 1L & d$area <= 11L))
  # every proportion within 4 standard errors of its rule
  near <- function(x, p) {
    expect_lt(abs(mean(x) - p), 4 * sqrt(p * (1 - p) / length(x)))
  }
  near(d$group == "r", 0.5)
  near(d$truth[d$group == "s"], 0.4)
  near(d$truth[d$group == "r"], 0.4)
  near(d$det[d$group == "r"], 0.4)
  near(d$det[d$group == "t"], 0.4 * 0.7)
  # the risk score: mean 250 + 100 bias off the reference, sd 40
  m <- tapply(d$risk, d$group, mean)
  expect_lt(abs(m[["r"]] - 250), 4 * 40 / sqrt(sum(d$group == "r")))
  expect_lt(abs(m[["t"]] - 280), 4 * 40 / sqrt(sum(d$group == "t")))
  expect_lt(abs(sd(d$risk[d$group == "s"]) - 40), 1.5)
  # reproducible from the seed; different seeds differ
  expect_identical(fairsim(500, seed = 2), fairsim(500, seed = 2))
  expect_false(identical(fairsim(500, seed = 2), fairsim(500, seed = 3)))
  set.seed(8)
  a <- fairsim(100)
  set.seed(8)
  expect_identical(fairsim(100), a)
  # clamping: a risk near the ends stays within 0 and 500, and a recording
  # rate above 1 is 1
  e <- fairsim(2000, base = 0.9, bias = -1, seed = 1)
  expect_true(all(e$det[e$group == "B"] == 1L))
  expect_true(all(e$risk >= 0 & e$risk <= 500))
  # a zero share: no record of that group
  expect_false(any(fairsim(300, props = c(1, 0), seed = 1)$group == "B"))
})

test_that("the audits recover the injected disparity", {
  # on seeded data the disparity injected lies within each audit's 99.9%
  # interval: the ratio of recorded rates 1 - bias, the score gap 100 bias
  for (s in 1:5) {
    b <- c(0.2, 0.5, -0.3, 0.6, 0)[s]
    d <- fairsim(6000, bias = b, seed = s)
    t <- fairdf(d$det, d$group, base = "A", conf = 0.999)$tab
    expect_lte(t$dilo[2], 1 - b)
    expect_gte(t$dihi[2], 1 - b)
    g <- scrdsp(d$risk, d$group, ref = "A", conf = 0.999)$tab
    expect_lte(g$lo[2], 100 * b)
    expect_gte(g$hi[2], 100 * b)
    # the true outcome does not depend on the group
    u <- fairdf(d$truth, d$group, base = "A", conf = 0.999)$tab
    expect_lte(u$dplo[2], 0)
    expect_gte(u$dphi[2], 0)
  }
  # the calibration audit sees the other group's areas ranked riskier than
  # their recorded outcomes, and the permutation test says so
  d <- fairsim(8000, nareas = 30, bias = 0.5, seed = 9)
  r <- calaud(d$area, d$risk, d$det, d$group, nperm = 999, seed = 1)
  expect_gt(r$groups$gap[2], 0)
  expect_lt(r$groups$p[2], 0.01)
  expect_lt(r$spear[["rho"]], 0)
  # with no disparity, no rank gap the test finds
  d <- fairsim(8000, nareas = 30, bias = 0, seed = 9)
  expect_gt(calaud(d$area, d$risk, d$det, d$group, nperm = 999,
                   seed = 1)$between[["p"]], 0.01)
})

test_that("nordet() counts officers in range and combines them", {
  r <- nordet(ref$nor$cx, ref$nor$ox, radius = 0.5, p = 0.7)
  expect_identical(r$k, as.integer(ref$nor$k))
  expect_equal(r$prob, ref$nor$prob, tolerance = 1e-15)
  set.seed(3)
  x <- data.frame(x = runif(50), y = runif(50))
  o <- cbind(runif(8), runif(8))
  r <- nordet(x, o, radius = 0.25, p = 0.4)
  dd <- sqrt(outer(x$x, o[, 1], "-")^2 + outer(x$y, o[, 2], "-")^2)
  k <- rowSums(dd <= 0.25)
  expect_identical(r$k, as.integer(k))
  expect_equal(r$prob, 1 - 0.6^k, tolerance = 1e-15)
  # every officer certain, no officers at all, and one exactly at range
  expect_identical(nordet(x, o, 0.25, p = 1)$prob, as.double(k > 0))
  expect_identical(nordet(x, NULL, 1)$k, integer(50))
  expect_identical(nordet(x, o[0, , drop = FALSE], 1)$prob, double(50))
  expect_identical(nordet(cbind(0, 0), cbind(3, 4), 5)$k, 1L)
  # draws: reproducible, near their probabilities
  d <- nordet(x, o, 0.25, 0.4, draw = TRUE, seed = 1)
  expect_identical(d, nordet(x, o, 0.25, 0.4, draw = TRUE, seed = 1))
  expect_named(d, c("k", "prob", "det"))
  expect_true(all(d$det[d$prob == 0] == 0L))
  big <- nordet(cbind(rep(0, 20000), 0), cbind(c(0, 0.1), 0), 1, 0.3,
                draw = TRUE, seed = 2)
  expect_lt(abs(mean(big$det) - 0.51), 4 * sqrt(0.51 * 0.49 / 20000))
  # a two-column table without the names x and y is read by position
  expect_identical(nordet(as.data.frame(o)[1:3, ], o, 0.2)$k,
                   nordet(o[1:3, ], o, 0.2)$k)
})

test_that("fairsim() and nordet() refuse what they cannot use", {
  errs("fairsim", list(
    list("whole number from 1", 0),
    list("two or more labels", groups = "A"),
    list("two or more labels", groups = c("A", NA)),
    list("not repeat a label", groups = c("A", "A")),
    list("must have length 2", props = c(1, 2, 3)),
    list("not all be 0", props = c(0, 0)),
    list("at least the number of groups", groups = 1:3, nareas = 2),
    list("from 0 to 1", base = 1.5),
    list("from -1 to 1", bias = 2)
  ))
  errs("nordet", list(
    list("two-column numeric matrix", 1:4, NULL, 1),
    list("columns x and y", data.frame(a = 1, b = 2, c = 3), NULL, 1),
    list("numeric coordinates", data.frame(x = "a", y = 1), NULL, 1),
    list("finite coordinates", cbind(NA, 1), NULL, 1),
    list("at least one point", matrix(0, 0, 2), NULL, 1),
    list("greater than 0", cbind(0, 0), NULL, 0),
    list("in \\(0, 1\\]", cbind(0, 0), NULL, 1, p = 0),
    list("TRUE or FALSE", cbind(0, 0), NULL, 1, draw = NA)
  ))
})
