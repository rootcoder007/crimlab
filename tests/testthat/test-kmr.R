# kmrec(): the curve, Greenwood's standard error, the interval, the median
# and the log-rank test against survival 3.8.6 (survfit(), survdiff(); the
# values below were computed once with it), hand cases, and recurrence.

gd <- data.frame(
  t = c(6, 6, 6, 7, 10, 13, 16, 22, 23, 6, 9, 10, 11, 17, 19, 20, 25, 32,
        32, 34, 35, 1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12,
        15, 17, 22, 23),
  e = c(rep(1, 9), rep(0, 12), rep(1, 21)),
  g = rep(c("drug", "placebo"), c(21, 21))
)

test_that("the curve and the test are survival's", {
  r <- kmrec(gd$t, gd$e, gd$g)
  k <- r$km[r$km$grp == "drug", ][1:8, ]
  expect_identical(k$time, c(6, 7, 9, 10, 11, 13, 16, 17))
  expect_identical(k$nrisk, c(21, 17, 16, 15, 13, 12, 11, 10))
  expect_equal(k$surv, c(0.85714285714285710, 0.80672268907563016,
                         0.80672268907563016, 0.75294117647058811,
                         0.75294117647058811, 0.69019607843137243,
                         0.62745098039215674, 0.62745098039215674),
               tolerance = 1e-14)
  expect_equal(k$se, c(0.076360354832121252, 0.086935285180057192,
                       0.086935285180057192, 0.096349652994320495,
                       0.096349652994320495, 0.106814707775009815,
                       0.114053865256752526, 0.114053865256752526),
               tolerance = 1e-13)
  expect_equal(k$lo, c(0.71981708391627042, 0.65312421846217106,
                       0.65312421846217106, 0.58591898202969417,
                       0.58591898202969417, 0.50961309910178021,
                       0.43939392496876695, 0.43939392496876695),
               tolerance = 1e-13)
  expect_equal(k$hi[1:2], c(1, 0.99644367590865912), tolerance = 1e-13)
  expect_identical(r$med$median, c(23, 8))
  expect_equal(r$lr$stat, 16.792940989216540, tolerance = 1e-12)
  expect_identical(r$lr$df, 1L)
  expect_equal(r$lr$p, pchisq(16.792940989216540, 1, lower.tail = FALSE),
               tolerance = 1e-10)
  g3 <- gd$g
  g3[c(1, 5, 30, 33)] <- "other"
  expect_equal(kmrec(gd$t, gd$e, g3)$lr$stat, 21.089742250036796,
               tolerance = 1e-12)
  expect_output(print(r), "log-rank")
})

test_that("hand cases: a median halfway, a curve that reaches 0", {
  # 1, 2, 3, 4: S = 0.75, 0.5, 0.25, 0; S sits at 1/2 from 2 to 3
  r <- kmrec(1:4)
  expect_equal(r$km$surv, c(0.75, 0.5, 0.25, 0))
  expect_identical(r$med$median, 2.5)
  expect_true(is.na(r$km$se[4]) && is.na(r$km$lo[4]))
  expect_null(r$lr)
  # at 1/2 with no later fall: the time it got there
  expect_identical(kmrec(c(1, 2, 3, 4), c(1, 1, 0, 0))$med$median, 2)
  expect_true(is.na(kmrec(c(1, 2, 3), c(1, 0, 0))$med$median))
  expect_output(print(r), "1 durations|4 durations")
})

test_that("recurrence: gaps within units, the last censored at end", {
  t <- c(2, 9, 30, 5, 70, 11, 12, 13)
  u <- rep(c("p", "q", "r"), c(3, 2, 3))
  r <- kmrec(t, unit = u, end = 100)
  ref <- kmrec(c(7, 21, 70, 65, 30, 1, 1, 87),
               event = c(1, 1, 0, 1, 0, 1, 1, 0))
  expect_identical(r$km, ref$km)
  expect_identical(r$n, 8L)
  # unordered times sort within each unit; no end, no censoring
  r2 <- kmrec(rev(t), unit = rev(u))
  expect_identical(r2$km$nevent, kmrec(c(7, 21, 65, 1, 1))$km$nevent)
  gg <- kmrec(t, grp = rep(c("x", "x", "x", "y", "y", "x", "x", "x")), unit = u)
  expect_identical(gg$med$grp, c("x", "y"))
})

test_that("bad input is refused", {
  errs("kmrec", list(
    list("`time` must be at least 0", c(1, -1)),
    list("`event`", 1:3, c(1, 2, 1)),
    list("not read with `unit`", 1:3, c(1, 1, 1), unit = c(1, 1, 1)),
    list("at or after every time", 1:3, unit = c(1, 1, 1), end = 2),
    list("same for every time", 1:3, grp = c(1, 2, 2), unit = c(1, 1, 2)),
    list("no gaps", 1:3, unit = 1:3),
    list("`level`", 1:3, level = 1)
  ))
})
