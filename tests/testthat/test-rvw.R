cpl <- function(k, n, level = 0.95) {
  as.numeric(stats::binom.test(k, n, conf.level = level)$conf.int)
}

test_that("reviews held and in time", {
  days <- c(12, 40, 75, 90, 130, 150, 64, 20, 61, 200)
  first <- c(NA, 31, 48, NA, 61, NA, 33, NA, 70, 59)
  r <- rvwtim(days, first, due = 60)
  expect_named(r, c("cov", "sum", "due", "level"))
  expect_identical(r$cov$band, c("up to 60", "over 60"))
  expect_identical(r$cov$n, c(3L, 7L))
  expect_identical(r$cov$never, c(2L, 2L))
  expect_equal(r$cov$share, c(2 / 3, 2 / 7), tolerance = 1e-15)
  expect_equal(c(r$cov$lo[2], r$cov$hi[2]), cpl(2, 7), tolerance = 1e-12)
  s <- r$sum
  expect_identical(s$group, "all")
  expect_identical(c(s$n, s$rev, s$past, s$pastnever, s$late),
                   c(10L, 6L, 7L, 2L, 2L))
  expect_equal(c(s$share, s$lo, s$hi), c(2 / 7, cpl(2, 7)),
               tolerance = 1e-12)
  expect_identical(s$med, stats::median(c(31, 48, 61, 33, 70, 59)))
  # groups in the order they appear, bands from cuts, another level
  g <- rvwtim(days, first, due = 60, by = rep(c("men", "women"), 5),
              cuts = c(30, 60, 120), level = 0.9)
  expect_identical(unique(g$cov$group), c("men", "women", "all"))
  expect_identical(unique(g$cov$band),
                   c("up to 30", "over 30 to 60", "over 60 to 120",
                     "over 120"))
  m <- g$cov[g$cov$group == "men", ]
  expect_identical(m$n, c(1L, 0L, 3L, 1L))
  expect_true(is.na(m$share[2]))
  expect_true(is.na(m$lo[2]))
  expect_true(is.na(m$hi[2]))
  expect_equal(c(m$lo[3], m$hi[3]), cpl(0, 3, 0.9), tolerance = 1e-12)
  # nobody reviewed: no median, and nothing past the day: no share
  z <- rvwtim(c(5, 8), c(NA, NA), due = 30)
  expect_true(is.na(z$sum$med))
  expect_true(is.na(z$sum$share))
  expect_identical(z$sum$past, 0L)
  expect_error(rvwtim(c(1, -2), c(1, 1), 5), "`days` must be numbers from 0")
  expect_error(rvwtim("a", 1, 5), "`days` must be numbers from 0")
  expect_error(rvwtim(c(1, 2), c(1, -1), 5), "`first` must be numbers")
  expect_error(rvwtim(c(1, 2), 1, 5), "`first` must be numbers")
  expect_error(rvwtim(c(1, 2), c("a", "b"), 5), "`first` must be numbers")
  expect_error(rvwtim(c(1, 2), c(1, 1), 0), "`due` must be one number")
  expect_error(rvwtim(c(1, 2), c(1, 1), 5, cuts = c(3, 3)),
               "`cuts` must be increasing numbers above 0\\.$")
  expect_error(rvwtim(c(1, 2), c(1, 1), 5, cuts = "a"), "`cuts` must be")
  expect_error(rvwtim(c(1, 2), c(1, 1), 5, by = "a"), "`by` must have one")
  expect_error(rvwtim(c(1, 2), c(1, 1), 5, level = 1), "`level` must be")
})

test_that("days on which an entitlement was missed", {
  miss <- c(100, 100, 80, 45, 10, 0, 100, 60, 20, 50)
  days <- c(20, 35, 16, 40, 18, 5, 62, 30, 17, 44)
  r <- entcmp(miss)
  expect_named(r, c("tab", "sum", "level"))
  expect_identical(r$tab$band, c("0", "over 0 to 20", "over 20 to 50",
                                 "over 50 to under 100", "100"))
  expect_identical(r$tab$n, c(1L, 2L, 2L, 2L, 3L))
  expect_equal(r$tab$share, c(1, 2, 2, 2, 3) / 10, tolerance = 1e-15)
  s <- r$sum
  expect_identical(c(s$n, s$none, s$all, s$half), c(10L, 1L, 3L, 6L))
  expect_equal(c(s$share, s$lo, s$hi), c(0.3, cpl(3, 10)),
               tolerance = 1e-12)
  expect_true(is.na(s$dayshr))
  f <- entcmp(miss, days = days, by = rep(c("east", "west"), each = 5),
              hours = c(0.4, 1.5, 2.2, 3, 3.5, NA, 0.8, 1, 0.5, 2),
              ref = c(10, 100, 5, 0, 0, 0, 30, 60, 51, 50),
              cuts = c(20, 50, 75), level = 0.99)
  expect_named(f, c("tab", "sum", "hrs", "ref", "level"))
  expect_identical(f$sum$group, c("east", "west", "all"))
  expect_equal(f$sum$dayshr[3], sum(miss * days) / (100 * sum(days)),
               tolerance = 1e-15)
  expect_equal(f$sum$dayshr[1],
               sum(miss[1:5] * days[1:5]) / (100 * sum(days[1:5])),
               tolerance = 1e-15)
  expect_equal(c(f$sum$lo[3], f$sum$hi[3]), cpl(3, 10, 0.99),
               tolerance = 1e-12)
  expect_identical(f$hrs$band, c("up to 0.5", "over 0.5 to 1", "over 1 to 2",
                                 "over 2 to 3", "over 3"))
  expect_identical(f$hrs$n, c(2L, 2L, 2L, 2L, 1L))
  expect_equal(sum(f$hrs$share), 1, tolerance = 1e-15)
  expect_identical(f$ref$most, c(1L, 2L, 3L))
  expect_identical(f$ref$net, c(4L, 1L, 5L))
  # one cut: one band between never and always
  o <- entcmp(c(0, 100, 33.3, 40), cuts = 33.3)
  expect_identical(o$tab$band, c("0", "over 0 to 33.3",
                                 "over 33.3 to under 100", "100"))
  expect_identical(o$tab$n, c(1L, 1L, 1L, 1L))
  expect_error(entcmp(c(10, 101)), "`miss` must be percentages")
  expect_error(entcmp(c(10, 20), days = c(1, 2, 3)),
               "`days` must have one value per stay")
  expect_error(entcmp(c(10, 20), days = c(1, -2)), "`days` must be numbers")
  expect_error(entcmp(c(10, 20), cuts = c(50, 100)),
               "`cuts` must be increasing numbers above 0 and below 100")
  expect_error(entcmp(c(10, 20), by = letters[1:3]), "`by` must have one")
  expect_error(entcmp(c(10, 20), level = 0), "`level` must be")
  expect_error(entcmp(c(10, 20), hours = c(1, -1)), "`hours` must be")
  expect_error(entcmp(c(10, 20), hours = 1), "`hours` must be")
  expect_error(entcmp(c(10, 20), ref = c(1, 200)), "`ref` must be perc")
  expect_error(entcmp(c(10, 20), ref = 5), "`ref` must have one value")
})

test_that("decision makers who do or do not decide alike", {
  dm <- rep(c("a", "b", "c", "d"), c(30, 28, 25, 4))
  dec <- rep(c("remain", "remove", "remain", "remove", "remain", "remove",
               "remain"), c(24, 6, 10, 18, 15, 10, 4))
  r <- dmvar(dec, dm, out = "remain")
  expect_named(r, c("tab", "test", "drop", "out", "level"))
  expect_identical(r$tab$dm, c("a", "b", "c", "d"))
  expect_identical(r$tab$n, c(30L, 28L, 25L, 4L))
  expect_identical(r$tab$k, c(24L, 10L, 15L, 4L))
  expect_equal(r$tab$rate, c(24 / 30, 10 / 28, 15 / 25, 1),
               tolerance = 1e-15)
  expect_equal(c(r$tab$lo[2], r$tab$hi[2]), cpl(10, 28), tolerance = 1e-12)
  expect_equal(sum(r$tab$share), 1, tolerance = 1e-15)
  ct <- suppressWarnings(stats::chisq.test(table(dm, dec), correct = FALSE))
  expect_equal(r$test$stat[1], unname(ct$statistic), tolerance = 1e-12)
  expect_identical(r$test$df[1], 3)
  expect_equal(r$test$p[1], ct$p.value, tolerance = 1e-12)
  expect_identical(r$test$low[1], as.numeric(sum(ct$expected < 5)))
  # without the decision maker of four cases
  k <- dm != "d"
  c3 <- stats::chisq.test(table(dm[k], dec[k]), correct = FALSE)
  expect_identical(r$test$set, c("all", "n >= nmin"))
  expect_identical(r$test$m, c(4L, 3L))
  expect_equal(r$test$stat[2], unname(c3$statistic), tolerance = 1e-12)
  expect_equal(r$drop$stat[4], unname(c3$statistic), tolerance = 1e-12)
  expect_equal(r$drop$p[4], c3$p.value, tolerance = 1e-12)
  # the other decision, a higher floor and another level
  v <- dmvar(dec, dm, out = "remove", nmin = 26, level = 0.9)
  expect_identical(v$out, "remove")
  expect_identical(v$test$m[2], 2L)
  expect_equal(v$test$stat[1], r$test$stat[1], tolerance = 1e-12)
  expect_equal(c(v$tab$lo[1], v$tab$hi[1]), cpl(6, 30, 0.9),
               tolerance = 1e-12)
  # two decision makers: none left to compare when one is dropped; and a
  # decision every one of them took
  w <- dmvar(c("y", "y", "y", "y"), c(1, 1, 2, 2), out = "y")
  expect_true(all(is.na(w$drop$stat)))
  expect_true(all(is.na(w$test$stat)))
  expect_identical(w$tab$dm, c("1", "2"))
  expect_error(dmvar(c("y", NA), c("a", "b"), "y"), "`dec` must be labels")
  expect_error(dmvar(c("y", "n"), "a", "y"), "`dm` must have one label")
  expect_error(dmvar(c("y", "n"), c("a", "b"), "maybe"), "`out` must be")
  expect_error(dmvar(c("y", "n"), c("a", "a"), "y"), "two or more decision")
  expect_error(dmvar(c("y", "n"), c("a", "b"), "y", nmin = 0.5),
               "`nmin` must be")
  expect_error(dmvar(c("y", "n"), c("a", "b"), "y", level = 2),
               "`level` must be")
})
