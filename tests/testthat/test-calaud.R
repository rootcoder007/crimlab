# calaud() and scrdsp(): against cor.test(), rank(), tapply(), an exact
# enumeration of the permutation test, oneway.test() and anova(lm()), the
# earlier package's recorded values, the disagreements with it, and
# refusals.

ref <- readRDS(test_path("ref-fair.rds"))
ca <- list(area = paste0("d", 1:8), risk = c(90, 80, 70, 60, 30, 20, 10, 5),
           outcome = c(10, 25, 20, 80, 70, 90, 85, 60),
           group = c("X", "X", "Y", "X", "Y", "Y", "X", "Y"))

test_that("the ranks, gaps and Spearman test follow base R", {
  set.seed(1)
  for (n in c(3, 5, 8, 9, 10, 14, 40)) {
    x <- rnorm(n)
    y <- x + rnorm(n)
    g <- rep(c("a", "b"), length.out = n)
    r <- calaud(seq_len(n), x, y, g, nperm = 0)
    o <- order(as.character(seq_len(n)))
    expect_identical(r$areas$riskrk, rank(-x)[o])
    expect_identical(r$areas$raterk, rank(-y)[o])
    expect_identical(r$areas$gap, (rank(-y) - rank(-x))[o])
    ct <- cor.test(x, y, method = "spearman")
    expect_equal(r$spear[["rho"]], unname(ct$estimate), tolerance = 1e-14)
    expect_equal(r$spear[["s"]], unname(ct$statistic), tolerance = 1e-12)
    # exact to nine areas, the Edgeworth series beyond: AS 89, as cor.test
    expect_equal(r$spear[["pex"]], ct$p.value, tolerance = 1e-12)
    cs <- cor.test(x, y, method = "spearman", exact = FALSE)
    expect_equal(r$spear[["pt"]], cs$p.value, tolerance = 1e-12)
    expect_false(r$ties)
    # the group means of the gaps, and the gaps sum to 0
    expect_equal(r$groups$gap, as.vector(tapply(rank(-y) - rank(-x), g, mean)),
                 tolerance = 1e-15)
    expect_equal(sum(r$areas$gap), 0, tolerance = 0)
  }
  # ties: average ranks, no exact p-value, the t one as cor.test gives it
  x <- c(1, 2, 2, 3, 4, 4, 5)
  y <- c(2, 1, 3, 3, 5, 6, 4)
  r <- calaud(letters[1:7], x, y, rep("g", 7), nperm = 0)
  expect_true(r$ties)
  expect_true(is.na(r$spear[["pex"]]))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  expect_equal(r$spear[["rho"]], unname(ct$estimate), tolerance = 1e-14)
  expect_equal(r$spear[["pt"]], ct$p.value, tolerance = 1e-12)
  expect_identical(r$areas$riskrk, rank(-x))
  # a constant ranking has no correlation
  r <- calaud(1:4, c(1, 1, 1, 1), 1:4, c(1, 1, 2, 2), nperm = 0)
  expect_true(all(is.nan(r$spear[c("rho", "pex", "t", "pt")])))
  # perfect agreement and disagreement
  expect_identical(calaud(1:5, 1:5, 1:5, rep(1, 5), nperm = 0)$spear[["pt"]],
                   0)
  expect_identical(calaud(1:5, 1:5, 5:1, rep(1, 5), nperm = 0)$spear[["rho"]],
                   -1)
})

test_that("records are rolled up to areas", {
  a <- c("a", "a", "b", "b", "c", "c", "c")
  rk <- c(10, 20, 30, 40, 5, 6, 7)
  oc <- c(1, 0, 1, 1, 0, 2, 1)
  gr <- c("x", "y", "y", "y", "x", "x", "y")
  pp <- c(a = 100, b = 200, c = 50)
  r <- calaud(a, rk, oc, gr, nperm = 0)
  expect_identical(r$areas$n, c(2L, 2L, 3L))
  expect_equal(r$areas$risk, as.vector(tapply(rk, a, mean)), tolerance = 0)
  expect_equal(r$areas$rate, as.vector(tapply(oc, a, mean)), tolerance = 0)
  # the most frequent group; a tie (area a: x and y once each) goes to the
  # first label, as table() and which.max() pick it
  expect_identical(r$areas$group, c("x", "y", "x"))
  expect_identical(r$areas$group, ref$agg$group)
  expect_identical(r$areas$risk, ref$agg$mean_risk)
  rp <- calaud(a, rk, oc, gr, population = pp, nperm = 0)
  expect_equal(rp$areas$rate, unname(as.vector(tapply(oc, a, sum)) / pp),
               tolerance = 1e-15)
  # the earlier package gives rates per 10,000 people
  expect_equal(rp$areas$rate * 1e4, ref$agg$outcome_rate, tolerance = 1e-14)
  # per record, or per area in any order, the same
  expect_identical(calaud(a, rk, oc, gr, pp[a], nperm = 0), rp)
  expect_identical(calaud(a, rk, oc, gr, rev(pp), nperm = 0), rp)
  # a data frame of the canonical columns
  df <- data.frame(area = a, risk = rk, outcome = oc, group = gr)
  expect_identical(calaud(df, nperm = 0), r)
  df$population <- pp[a]
  expect_identical(calaud(df, nperm = 0), rp)
  # named vectors are matched by name
  nm <- paste0("r", 1:7)
  o <- c(3:7, 1:2)
  expect_identical(
    calaud(setNames(a, nm), setNames(rk, nm)[o], setNames(oc, nm),
           setNames(gr, nm)[rev(o)], nperm = 0), r)
})

test_that("the values agree with the earlier package's", {
  r <- calaud(ca$area, ca$risk, ca$outcome, ca$group, nperm = 0)
  expect_equal(r$spear[["rho"]], ref$cal$spearman, tolerance = 1e-14)
  expect_equal(r$spear[["pex"]], ref$cal$spearman_pvalue, tolerance = 1e-14)
  expect_equal(r$groups$gap, unname(unlist(ref$cal$group_rank_gap)[
    r$groups$group]), tolerance = 1e-15)
  expect_equal(r$areas$gap, unname(unlist(ref$cal$rank_gap)[r$areas$area]),
               tolerance = 0)
  # the earlier package's other copy of the audit uses the large-sample
  # p-value only; calaud() gives both, as pex and pt
  expect_equal(r$spear[["pt"]], ref$cal2$spearman_pvalue, tolerance = 1e-14)
})

test_that("the permutation p-values match an exact enumeration", {
  r <- calaud(ca$area, ca$risk, ca$outcome, ca$group, nperm = 20000, seed = 1)
  gap <- r$areas$gap
  obs <- sum(gap[r$areas$group == "X"])
  ss <- function(s) s^2 / 4 + (-s)^2 / 4
  # every way of choosing the 4 areas of group X among the 8
  s <- combn(8, 4, function(i) sum(gap[i]))
  pe <- mean(abs(s) >= abs(obs) - 1e-9)
  pss <- mean(ss(s) >= ss(obs) - 1e-9)
  # a Monte Carlo p-value of B draws is within 4 standard errors of the
  # exact one (the +1 in the numerator and denominator moves it by < 1/B)
  se <- sqrt(pe * (1 - pe) / 20000)
  expect_lt(abs(r$groups$p[1] - pe), 4 * se + 1e-4)
  expect_lt(abs(r$between[["p"]] - pss), 4 * se + 1e-4)
  expect_equal(r$between[["ss"]], ss(obs), tolerance = 1e-14)
  # reproducible from the seed, and from set.seed()
  cd <- as.data.frame(ca)
  expect_identical(calaud(cd, nperm = 99, seed = 2),
                   calaud(ca$area, ca$risk, ca$outcome, ca$group, nperm = 99,
                          seed = 2))
  set.seed(4)
  a1 <- calaud(cd, nperm = 99)
  set.seed(4)
  expect_identical(calaud(cd, nperm = 99), a1)
  expect_gte(min(a1$groups$p), 1 / 100)
  # one group holding every area: its mean gap is 0, its p-value 1
  r <- calaud(1:5, c(5, 3, 4, 1, 2), 1:5, rep("g", 5), nperm = 50, seed = 1)
  expect_identical(r$groups$p, 1)
  # a label with no area of its own (every area's majority is another)
  r <- calaud(c(1, 1, 2, 2, 3, 3), 1:6, 6:1, c("x", "x", "x", "x", "x", "y"),
              nperm = 10, seed = 1)
  expect_identical(r$groups$n, c(3L, 0L))
  expect_true(is.nan(r$groups$gap[2]))
  expect_true(is.nan(r$groups$p[2]))
})

test_that("scrdsp() is the one-way analysis of variance", {
  set.seed(6)
  g <- rep(c("a", "b", "c"), c(12, 9, 15))
  x <- rnorm(36, c(a = 0, b = 0.4, c = 1)[g])
  r <- scrdsp(x, g, ref = "a", conf = 0.9)
  t <- r$tab
  expect_equal(t$mean, as.vector(tapply(x, g, mean)), tolerance = 1e-15)
  expect_equal(t$median, as.vector(tapply(x, g, median)), tolerance = 1e-15)
  expect_equal(t$sd, as.vector(tapply(x, g, sd)), tolerance = 1e-14)
  ow <- oneway.test(x ~ g, var.equal = TRUE)
  expect_equal(unname(r$anova), unname(c(ow$statistic, ow$parameter,
                                         ow$p.value)), tolerance = 1e-12)
  wl <- oneway.test(x ~ g)
  expect_equal(unname(r$welch), unname(c(wl$statistic, wl$parameter,
                                         wl$p.value)), tolerance = 1e-12)
  # gaps with pooled-variance intervals
  expect_identical(r$ref, "a")
  fit <- lm(x ~ g)
  ci <- confint(fit, level = 0.9)
  expect_equal(t$gap[2:3], unname(coef(fit)[2:3]), tolerance = 1e-13)
  expect_equal(t$se[2:3], unname(sqrt(diag(vcov(fit)))[2:3]),
               tolerance = 1e-13)
  expect_equal(cbind(t$lo, t$hi)[2:3, ], unname(ci[2:3, ]), tolerance = 1e-13)
  expect_equal(r$spread, diff(range(t$mean)), tolerance = 1e-15)
  expect_identical(scrdsp(x, g, ref = "b")$tab$gap[2], 0)
  # by default the reference is the group with the lowest mean
  expect_identical(scrdsp(c(3, 4, 1, 2), c("p", "p", "q", "q"))$ref, "q")
})

test_that("scrdsp() with a group of one, and the earlier package's", {
  x <- c(9, 10, 11, 19, 20, 21, 30, 14)
  g <- c("A", "A", "A", "B", "B", "B", "C", "A")
  r <- scrdsp(x, g)
  # one case: no standard deviation, no Welch test; the analysis of
  # variance is anova(lm()), to which the case adds between-group spread
  expect_true(is.nan(r$tab$sd[3]))
  expect_true(all(is.nan(r$welch)))
  a <- anova(lm(x ~ g))
  expect_equal(r$anova[["f"]], a$`F value`[1], tolerance = 1e-13)
  expect_equal(r$anova[["p"]], a$`Pr(>F)`[1], tolerance = 1e-12)
  expect_equal(r$tab$mean, unname(unlist(ref$sd$group_means)),
               tolerance = 1e-15)
  expect_equal(r$tab$gap, unname(unlist(ref$sd$gaps)), tolerance = 1e-15)
  expect_identical(r$spread, ref$sd$spread)
  expect_identical(r$ref, ref$sd$reference)
  # the earlier package gives no F here (oneway.test() refuses a group of
  # one) in one copy, and drops the group in the other
  expect_true(is.na(ref$sd$anova_f))
  k <- g != "C"
  d2 <- scrdsp(x[k], g[k])$anova
  expect_equal(d2[["f"]], ref$sd2$anova_f, tolerance = 1e-13)
  expect_equal(d2[["p"]], ref$sd2$anova_pvalue, tolerance = 1e-12)
  # every case its own group: no within-group variation to compare with
  r <- scrdsp(1:3, c("a", "b", "c"))
  expect_true(is.nan(r$anova[["f"]]))
  expect_true(is.nan(r$tab$lo[2]))
})

test_that("calaud() and scrdsp() refuse what they cannot use", {
  errs("calaud", list(
    list("are needed", 1:3, 1:3),
    list("not both", data.frame(area = 1:3), risk = 1:3),
    list("needs columns", data.frame(area = 1:3, risk = 1:3)),
    list("at least one record", character(0), numeric(0), numeric(0),
         character(0)),
    list("at least 3 areas", c(1, 2, 2), 1:3, 1:3, 1:3),
    list("finite numbers", 1:3, c(1, NA, 3), 1:3, 1:3),
    list("must have length 3", 1:3, 1:3, 1:2, 1:3),
    list("greater than 0", 1:3, 1:3, 1:3, 1:3, population = c(1, 0, 1)),
    list("one value within each area", c(1, 1, 2, 3), 1:4, 1:4, 1:4,
         population = c(1, 2, 3, 4)),
    list("no value for area 3", 1:3, 1:3, 1:3, 1:3,
         population = c(`1` = 1, `2` = 2, `4` = 3)),
    list("names an area twice", 1:3, 1:3, 1:3, 1:3,
         population = setNames(1:3, c(1, 1, 3))),
    list("one per area named", 1:3, 1:3, 1:3, 1:3, population = c(1, 2)),
    list("whole number", 1:3, 1:3, 1:3, 1:3, nperm = 1.5)
  ))
  errs("scrdsp", list(
    list("at least two groups", 1:3, c(1, 1, 1)),
    list("finite numbers", c(1, Inf), 1:2),
    list("names no group", 1:4, c(1, 1, 2, 2), ref = "z"),
    list("strictly between", 1:4, c(1, 1, 2, 2), conf = 0)
  ))
})
