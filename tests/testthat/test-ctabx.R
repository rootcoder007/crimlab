# ctabx(): the statistics against stats::chisq.test() and fisher.test(),
# Cramer's V and its interval against DescTools, the measures by hand,
# aggregation, empty and degenerate tables, the input forms and the
# refusals.

set.seed(7)
ct <- data.frame(
  a = sample(c("x", "y", "z"), 120, TRUE),
  b = sample(c("p", "q", "r", "s"), 120, TRUE)
)

test_that("Pearson and likelihood-ratio statistics are chisq.test()'s", {
  r <- ctabx(ct$a, ct$b)
  ref <- suppressWarnings(chisq.test(table(ct$a, ct$b), correct = FALSE))
  s <- r$stat
  expect_equal(s$chisq, unname(ref$statistic), tolerance = 1e-12)
  expect_identical(s$df, 6)
  expect_equal(s$p, ref$p.value, tolerance = 1e-10)
  expect_equal(unname(r$expd), unname(ref$expected), tolerance = 1e-13)
  expect_equal(unname(r$res), unname(ref$residuals), tolerance = 1e-12)
  o <- table(ct$a, ct$b)
  e <- ref$expected
  expect_equal(s$g2, 2 * sum(o * log(o / e)), tolerance = 1e-12)
  expect_equal(s$v, sqrt(s$chisq / (120 * 2)), tolerance = 1e-14)
  expect_equal(s$cc, sqrt(s$chisq / (s$chisq + 120)), tolerance = 1e-14)
  expect_identical(
    dimnames(r$tab),
    list(c("x", "y", "z"), c("p", "q", "r", "s"))
  )
})

test_that("Cramer's V and its interval are DescTools'", {
  skip_if_not_installed("DescTools")
  r <- ctabx(ct$a, ct$b)$stat
  d <- DescTools::CramerV(table(ct$a, ct$b), conf.level = 0.95)
  # DescTools stops its root search when the probability is within 1e-5
  expect_equal(c(r$v, r$vlo, r$vhi), unname(d), tolerance = 1e-3)
  b <- DescTools::CramerV(table(ct$a, ct$b), correct = TRUE)
  expect_equal(r$vbc, b, tolerance = 1e-12)
})

test_that("2 x 2: Yates, phi and Fisher's test are R's", {
  m <- matrix(c(12, 5, 9, 21), 2, dimnames = list(c("a", "b"), c("y", "n")))
  r <- ctabx(m, correct = TRUE, exact = TRUE)
  s <- r$stat
  ref <- chisq.test(m)
  expect_equal(s$chisq, unname(ref$statistic), tolerance = 1e-12)
  f <- fisher.test(m)
  expect_equal(s$fp, f$p.value, tolerance = 1e-10)
  # fisher.test() finds its roots in t on (0, 1) to uniroot()'s 1.2e-4,
  # and the upper end is 1 / t with t near 0.04 here: up to 3e-3 relative
  expect_equal(
    c(s$or, s$orlo),
    unname(c(f$estimate, f$conf.int[1])),
    tolerance = 1e-4
  )
  expect_equal(s$orhi, f$conf.int[2], tolerance = 5e-3)
  expect_equal(
    s$phi,
    (12 * 21 - 9 * 5) / sqrt(21 * 26 * 17 * 30),
    tolerance = 1e-14
  )
  # V is always from the uncorrected statistic
  expect_equal(s$v, abs(s$phi), tolerance = 1e-14)
  # the odds ratio inverts when the rows swap
  r2 <- ctabx(m[2:1, ], exact = TRUE)$stat
  expect_equal(r2$or, 1 / s$or, tolerance = 1e-9)
  expect_equal(r2$orlo, 1 / s$orhi, tolerance = 1e-9)
  # a zero cell: the estimate at the edge of its range
  z <- ctabx(matrix(c(0, 5, 9, 21), 2), exact = TRUE)$stat
  expect_identical(c(z$or, z$orlo), c(0, 0))
  z2 <- ctabx(matrix(c(4, 0, 9, 21), 2), exact = TRUE)$stat
  expect_identical(z2$orhi, Inf)
})

test_that("r x c Fisher's test is fisher.test()'s", {
  m <- matrix(c(3, 1, 0, 2, 4, 1, 0, 2, 5), 3)
  r <- ctabx(m, exact = TRUE)$stat
  expect_equal(r$fp, fisher.test(m)$p.value, tolerance = 1e-10)
  expect_true(is.na(r$or))
  expect_error(
    ctabx(matrix(c(80, 75, 90, 60, 70, 85, 77, 66, 91), 3), exact = TRUE),
    "more than 1e7 tables"
  )
})

test_that("aggregation, distinct ids, the diagonal and empty rows", {
  d <- data.frame(
    r = c("A", "A", "B", "B", "B"),
    c = c("F", "M", "F", "F", "M"),
    n = c(3, 4, 5, 6, 7),
    id = c(1, 2, 1, 1, 3)
  )
  expect_identical(unname(ctabx(d$r, d$c, d$n)$tab), matrix(c(3, 11, 4, 7), 2))
  expect_identical(
    unname(ctabx(d$r, d$c, d$n, "max")$tab),
    matrix(c(3, 6, 4, 7), 2)
  )
  expect_identical(
    unname(ctabx(d$r, d$c, d$n, "min")$tab),
    matrix(c(3, 5, 4, 7), 2)
  )
  expect_identical(
    unname(ctabx(d$r, d$c, d$n, "mean")$tab),
    matrix(c(3, 5.5, 4, 7), 2)
  )
  expect_identical(
    unname(ctabx(d$r, d$c, d$id, "dist")$tab),
    matrix(c(1, 1, 1, 1), 2)
  )
  expect_true(is.na(ctabx(d$r, d$c, d$n, "max")$stat$chisq))
  # the diagonal: labels matched by name, not position
  g <- ctabx(c("N", "N", "E", "W", "E"), c("N", "E", "E", "W", "E"))
  expect_identical(g$stat$diag, 4 / 5)
  expect_true(is.na(ctabx(ct$a, ct$b)$stat$diag))
  # an unused level is dropped from the statistics, kept in the table
  f <- factor(c("u", "u", "v", "v"), levels = c("u", "v", "w"))
  h <- ctabx(f, c("a", "b", "a", "a"))
  expect_identical(h$drop, 0L)
  m <- matrix(c(5, 0, 3, 3, 0, 4), 3)
  k <- ctabx(m)
  expect_identical(k$drop, 1L)
  expect_identical(dim(k$expd), c(2L, 2L))
  # one row: no test
  one <- ctabx(c("a", "a"), c("x", "y"))
  expect_true(is.na(one$stat$chisq))
  expect_output(print(ctabx(ct$a, ct$b)), "Cramer's V")
  expect_output(print(k), "empty rows or columns")
  expect_output(print(ctabx(m[c(1, 3), ], exact = TRUE)), "Fisher")
  expect_output(print(g), "diagonal share")
})

test_that("named vectors match by name; a clpre list fills the arguments", {
  a <- stats::setNames(ct$a, seq_len(120))
  b <- stats::setNames(ct$b, seq_len(120))
  expect_identical(ctabx(a[120:1], b)$tab, ctabx(ct$a, ct$b)$tab)
  p <- structure(list(row = ct$a, col = ct$b), class = "clpre")
  expect_identical(ctabx(p)$stat, ctabx(ct$a, ct$b)$stat)
})

test_that("refusals", {
  errs(
    "ctabx",
    list(
      list("`row` must be the row labels", 1:3),
      list("`w` goes with", matrix(1, 2, 2), NULL, 1),
      list("`row` must hold counts", matrix(c(1, -1, 2, 2), 2)),
      list("`agg` applies to", matrix(1, 2, 2), agg = "max"),
      list("leave `w` out", c("a", "b"), c("x", "y"), c(1, 2), "n"),
      list("needs `w`", c("a", "b"), c("x", "y"), agg = "sum"),
      list("`w` must be at least 0", c("a", "b"), c("x", "y"), c(1, -2)),
      list("`col` must be labels", c("a", "b"), c("x", NA)),
      list("whole counts", matrix(c(1.5, 2, 3, 4), 2), exact = TRUE),
      list(
        "needs counts",
        c("a", "b"),
        c("x", "y"),
        c(1, 2),
        "max",
        exact = TRUE
      ),
      list("`agg` must be one of", c("a", "b"), c("x", "y"), agg = "median"),
      list("`exact` must be TRUE or FALSE", matrix(1, 2, 2), exact = NA),
      list("`row` must have at least one label", character(0), character(0))
    )
  )
})

# Parity with rmorie (helper-port.R): live when installed, recorded otherwise.

test_that("ctabx() reproduces the OTIS cross-tabulations", {
  b <- portref("otisxtab")
  chk <- function(r, ref) {
    for (row in ref$rows) {
      expect_identical(
        as.character(r$tab[row[1], ref$head[2:(length(row) - 1)]]),
        row[2:(length(row) - 1)]
      )
    }
  }
  c4 <- pin$c04
  chk(
    ctabx(
      c4$Race,
      c4$Region_MostRecentPlacement,
      c4$NumberIndividuals_RestrictiveConfinement
    ),
    b$c04
  )
  b4 <- pin$b04
  chk(
    ctabx(
      b4$Region_AtTimeOfPlacement,
      b4$Measure,
      b4$NumberConsecutiveDays_Segregation,
      "max"
    ),
    b$b04
  )
  d2 <- pin$d02
  yr <- ctabx(d2$Year, rep("all", nrow(d2)), d2$Number_CustodialDeaths)$tab
  expect_identical(
    as.character(yr[, 1]),
    vapply(b$d02[[1]], function(e) e[2], "")
  )
  g <- ctabx(d2$Gender, d2$Year, d2$Number_CustodialDeaths)$tab
  expect_identical(as.character(g["Male", ]), b$d02[[2]][[2]][2:7])
})

test_that("ctabx() reproduces the c-series chi-square family", {
  b <- portref("cchi2")
  cs <- pin$cser
  sl <- list(
    c03 = c("Race", "Gender"),
    c04 = c("Race", "Region_MostRecentPlacement"),
    c05 = c("Religion", "Region_MostRecentPlacement"),
    c06 = c("Age_Category", "Region_MostRecentPlacement"),
    c07 = c("Alert_Type", "Gender"),
    c08 = c("Religion", "Gender"),
    c09 = c("Age_Category", "Gender")
  )
  for (k in names(sl)) {
    s <- ctabx(
      cs[[sl[[k]][1]]],
      cs[[sl[[k]][2]]],
      cs$NumberIndividuals_RestrictiveConfinement,
      correct = TRUE
    )$stat
    expect_equal(
      c(s$chisq, s$df, s$p),
      unname(b[[k]][1:3]),
      tolerance = 1e-12,
      info = k
    )
    if (k != "c03") {
      expect_equal(s$v, unname(b[[k]][4]), tolerance = 1e-12, info = k)
    }
  }
  # disagreement: on a 2 x 2 table rmorie takes V from the continuity-
  # corrected statistic; crimlab's V is always from Pearson's
  s <- ctabx(
    cs$Race,
    cs$Gender,
    cs$NumberIndividuals_RestrictiveConfinement,
    correct = TRUE
  )$stat
  expect_equal(unname(b$c03[4]), sqrt(s$chisq / s$n), tolerance = 1e-12)
  expect_gt(s$v, unname(b$c03[4]))
})

test_that("ctabx() reproduces locality, weapons and co-occurrence", {
  b <- portref("locality")
  l <- pin$loc
  s <- ctabx(l$Region_AtTimeOfPlacement, l$Region_MostRecentPlacement)$stat
  # mrm_otis_region_locality() rounds what it returns
  expect_identical(
    c(round(s$chisq, 2), s$df, signif(s$p, 3), round(s$v, 4), round(s$diag, 4)),
    b$otis
  )
  expect_equal(c(s$chisq, s$p, s$v, s$diag), b$uof5, tolerance = 1e-12)
  # disagreement: on a 2 x 2 table mrm_uof_region_locality() takes V from
  # the continuity-corrected statistic
  l2 <- pin$loc2
  s2 <- ctabx(l2$at, l2$now, correct = TRUE)$stat
  expect_equal(
    c(s2$chisq, s2$p, s2$diag),
    b$uof2[c(1, 2, 4)],
    tolerance = 1e-12
  )
  expect_equal(b$uof2[3], sqrt(s2$chisq / s2$n), tolerance = 1e-12)
  w <- portref("weapon")
  r <- ctabx(pin$wpn$weapon, pin$wpn$force)
  expect_equal(
    c(r$stat$chisq, r$stat$df, r$stat$p, r$stat$v),
    w$stat,
    tolerance = 1e-12
  )
  expect_equal(
    sort(abs(r$res), decreasing = TRUE)[1:3],
    abs(w$top),
    tolerance = 1e-12
  )
  m <- portref("mort")
  mo <- pin$mort
  for (i in 1:3) {
    s <- ctabx(mo[[m$alert_a[i]]] == "Yes", mo[[m$alert_b[i]]] == "Yes")$stat
    expect_identical(
      c(round(s$chisq, 2), signif(s$p, 3), round(s$v, 4)),
      c(m$chi2[i], m$p_value[i], m$cramers_v[i])
    )
  }
})

test_that("ctabx() reproduces the SIU tables", {
  b <- portref("siuchi")
  y <- pin$siuy
  yr <- substr(y$date_of_incident_iso, 1, 4)
  s <- ctabx(yr, y$charges_recommended, correct = TRUE)
  expect_equal(
    c(s$stat$chisq, s$stat$df, s$stat$p),
    unname(b$stat),
    tolerance = 1e-12
  )
  expect_equal(rowSums(s$tab), b$years, tolerance = 0)
  o <- ctabx(pin$siuo$police_service, pin$siuo$director_decision_category)$tab
  ref <- b$out
  got <- o[cbind(as.character(ref$service), as.character(ref$outcome))]
  expect_equal(got, ref$n_cases, tolerance = 0)
  expect_equal(
    round(
      (o / rowSums(o))[cbind(
        as.character(ref$service),
        as.character(ref$outcome)
      )],
      4
    ),
    ref$share_within_service,
    tolerance = 0
  )
})
