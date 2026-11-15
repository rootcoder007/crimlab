# Statistics Canada's published tables, the gaps against a reference, the
# generations of exposure and the Calls to Action: the source package's
# tests of these functions, then the values that package gives.

test_that(
  paste0(
    "sctidy() keeps the quality flags and ",
    "never keeps a suppressed value"
  ),
  {
    d <- data.frame(
      GEO = "Canada",
      VALUE = c("45.8", "12", "", "3", "7"),
      STATUS = c("", "E", "x", "F", ".."),
      stringsAsFactors = FALSE
    )
    t <- sctidy(d)
    expect_s3_class(t, "sctidy")
    expect_identical(as.character(t$quality), c("ok", "E", "x", "F", ".."))
    expect_identical(t$VALUE, c(45.8, 12, NA, NA, NA))
    expect_error(
      sctidy(data.frame(a = 1)),
      "no VALUE column, nor census columns"
    )
    expect_error(sctidy(1), "must be a data frame")
  }
)


test_that("sctidy() makes a census table's wide layout long", {
  w <- data.frame(
    REF_DATE = 2021,
    GEO = "Canada",
    `Age (20)` = c("Total - Age", "0 to 14 years"),
    Coordinate = c("1.1", "1.2"),
    `Indigenous identity (9):First Nations (North American Indian)[4]` = c(
      1048405,
      300000
    ),
    Symbol = c("", "x"),
    `Indigenous identity (9):Inuk (Inuit)[6]` = c(70540, 22000),
    Symbol2 = c("E", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # StatCan repeats the Symbol column after every value column
  names(w)[names(w) == "Symbol2"] <- "Symbol"
  t <- sctidy(w)
  expect_identical(
    t[["Indigenous identity (9)"]],
    rep(c("First Nations (North American Indian)", "Inuk (Inuit)"), each = 2)
  )
  expect_identical(t$VALUE, c(1048405, NA, 70540, 22000))
  expect_identical(as.character(t$quality), c("ok", "x", "E", "ok"))
  expect_identical(t[["Age (20)"]], rep(c("Total - Age", "0 to 14 years"), 2))
  e <- scests(t)
  expect_identical(nrow(e), 4L)
  expect_false("Coordinate" %in% names(e))
  # one label, two members: told apart by the member number in the coordinate
  du <- data.frame(
    REF_DATE = 2021,
    GEO = "Canada",
    DGUID = "2021A000011124",
    `Household type (16)` = c("Without children", "Without children", "Total"),
    Coordinate = c("1.3", "1.7", "1.1"),
    `Tenure (8):Owner[2]` = c(5, 6, 11),
    Symbol = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  tu <- sctidy(du)
  expect_identical(
    tu[["Household type (16)"]],
    c("Without children [3]", "Without children [7]", "Total")
  )
  expect_identical(nrow(scests(tu)), 3L)
  bad <- data.frame(`A (2):x[1]` = 1, `B (2):y[1]` = 2, check.names = FALSE)
  expect_error(sctidy(bad), "more than one dimension")
})


test_that(
  paste0(
    "indgap() works the ratio and ",
    "difference from the published intervals"
  ),
  {
    d <- data.frame(
      g = c(
        "First Nations",
        "Metis",
        "Total Indigenous",
        "Inuit",
        "Non-Indigenous"
      ),
      est = c(45.8, 52.8, 49, 47.8, 61.5),
      lo = c(43.9, 50.6, 47.5, 43.1, 61),
      hi = c(47.7, 55, 50.5, 52.5, 62),
      q = c("", "", "", "E", ""),
      stringsAsFactors = FALSE
    )
    r <- indgap(d, "g", "est", "lo", "hi", qual = "q")
    z <- stats::qnorm(0.975)
    se <- function(l, h) (h - l) / (2 * z)
    fn <- r[r$group == "First Nations", ]
    expect_identical(fn$difference, 45.8 - 61.5)
    expect_identical(
      fn$diflo,
      (45.8 - 61.5) - z * sqrt(se(43.9, 47.7)^2 + se(61, 62)^2)
    )
    sl <- sqrt(
      ((log(47.7) - log(43.9)) / (2 * z))^2 + ((log(62) - log(61)) / (2 * z))^2
    )
    expect_identical(fn$ratioup, exp(log(45.8 / 61.5) + z * sl))
    expect_true(r$usewithcau[r$group == "Inuit"])
    expect_identical(r$poogrp, c(FALSE, FALSE, TRUE, FALSE))
    d$q[2] <- "F"
    expect_true(is.na(indgap(d, "g", "est", "lo", "hi", qual = "q")$ratio[2]))
    expect_error(
      indgap(d, "g", "est", "lo", "hi", ref = "zz"),
      "exactly one row"
    )
    d2 <- d
    d2$lo[1] <- 50
    expect_warning(
      {
        r2 <- indgap(d2, "g", "est", "lo", "hi")
      },
      "does not contain the estimate for: First Nations"
    )
    expect_true(is.na(r2$difference[1]))
    d4 <- d
    d4$lo[5] <- 62.5
    expect_error(
      indgap(d4, "g", "est", "lo", "hi"),
      "reference group's published interval"
    )
    d3 <- d
    d3$q <- c("", "", "", "", "x")
    expect_error(
      indgap(d3, "g", "est", "lo", "hi", qual = "q"),
      "reference estimate is suppressed"
    )
    expect_error(indgap(d, "zz", "est", "lo", "hi"), "column `zz` not found")
    expect_error(
      indgap(d, "g", "est", "lo", "hi", qual = "zz"),
      "column `zz` not found"
    )
    expect_error(
      indgap(1, "g", "est", "lo", "hi"),
      "`data` must be a data frame"
    )
  }
)


test_that("genexpo() codes 0, 1, 2 and keeps what is not known unknown", {
  expect_identical(
    genexpo(
      c(TRUE, FALSE, TRUE, NA, FALSE, NA, TRUE),
      c(TRUE, FALSE, FALSE, TRUE, NA, NA, NA)
    ),
    c(2L, 0L, 1L, NA, NA, NA, NA)
  )
  expect_identical(
    genexpo(
      c(TRUE, FALSE, TRUE, NA, FALSE, NA, TRUE),
      c(TRUE, FALSE, FALSE, TRUE, NA, NA, NA),
      partial = "lowerbnd"
    ),
    c(2L, 0L, 1L, 1L, NA, NA, 1L)
  )
  expect_error(genexpo(1, TRUE), "`parent` must be TRUE, FALSE or NA")
  expect_error(genexpo(TRUE, c(TRUE, FALSE)), "same length")
})


test_that(
  paste0(
    "exptrnd(): odds ratios by hand, the trend test as prop.trend.test ",
    "gives it"
  ),
  {
    r <- exptrnd(c(60, 95, 70), c(500, 520, 300))
    expect_identical(r$orratio[3], (70 / 230) / (60 / 440))
    se <- sqrt(1 / 70 + 1 / 230 + 1 / 60 + 1 / 440)
    expect_identical(
      r$lower[3],
      exp(log((70 / 230) / (60 / 440)) - stats::qnorm(0.975) * se)
    )
    expect_equal(
      attr(r, "trend")$p.value,
      stats::prop.trend.test(c(60, 95, 70), c(500, 520, 300), 0:2)$p.value,
      tolerance = 1e-12
    )
    expect_no_warning({
      r0 <- exptrnd(c(0, 5), c(10, 10))
    })
    expect_identical(r0$orratio[2], (5.5 / 5.5) / (0.5 / 10.5))
    expect_error(exptrnd(c(1, 20), c(10, 10)), "0 <= cases <= n")
    expect_error(exptrnd(1, 10), "two levels at least")
    expect_error(
      exptrnd(c(1, 2), c(10, 10), dose = c(1, 0)),
      "strictly increasing"
    )
  }
)


test_that(
  paste0(
    "trccmpl() and trcprog(): the counts the reports state, the exact rate ",
    "interval"
  ),
  {
    tc <- trccmpl()
    expect_identical(tc$completed, c(8L, 0L, 3L, 2L))
    expect_identical(tail(tc$cumulative, 1), 13L)
    p <- trcprog()
    span <- 2023 - 2015.96
    expect_identical(p$rate, 13 / span)
    expect_equal(
      p$ratelo,
      stats::qgamma(0.025, 13) / span,
      tolerance = testthat::testthat_tolerance()
    )
    expect_equal(p$rateup, stats::qgamma(0.975, 14) / span, tolerance = 1e-14)
    expect_identical(p$cmplby, 2023 + 81 / (13 / span))
    expect_identical(trcprog(94, 2030)$cmplby, 2031)
    expect_identical(trcprog(0, 2020)$cmplbylate, Inf)
    expect_error(trcprog(c(5, 3), c(2020, 2021)), "cumulative counts")
    expect_error(trcprog(1, c(2020, 2021)), "one count per year")
    expect_error(trcprog(1, 2010), "after `start`")
  }
)


test_that("qltyrcd() refuses non-character flags and warns of unknown ones", {
  expect_error(qltyrcd("a", suppr = mean), "`suppr` must be data")
  expect_error(
    qltyrcd(c(1, 2)),
    "^qltyrcd: `x` must be a vector of flags \\(character\\)"
  )
  expect_error(qltyrcd(TRUE), "vector of flags \\(character\\)")
  fl <- c("np", "*", "Q", "", NA)
  expect_warning(
    qltyrcd(fl, suppr = "np", caution = "*"),
    "\"Q\" are in neither `suppr` nor `caution`"
  )
  r <- suppressWarnings(qltyrcd(fl, suppr = "np", caution = "*"))
  expect_identical(r, c("x", "E", "ok", "ok", "ok"))
  expect_silent(qltyrcd(factor(c("E", "F", "ok"))))
  expect_identical(qltyrcd(NA), "ok")
})


test_that("indgap() refuses a `pooled` label that is not a group", {
  d <- data.frame(
    g = c("A", "B", "Total", "R"),
    e = c(10, 12, 11, 8),
    lo = c(8, 10, 9, 7),
    hi = c(12, 14, 13, 9),
    stringsAsFactors = FALSE
  )
  expect_error(
    indgap(d, "g", "e", "lo", "hi", ref = "R", pooled = "Totl"),
    "^indgap: `pooled` names \"Totl\", not a group of `data`"
  )
  expect_error(
    indgap(d, "g", "e", "lo", "hi", ref = "R", pooled = mean),
    "`pooled` must be data"
  )
  r <- indgap(d, "g", "e", "lo", "hi", ref = "R", pooled = "Total")
  expect_identical(r$poogrp, c(FALSE, FALSE, TRUE))
})


test_that("refusals name the function", {
  expect_error(limfrom(1), "^limfrom: give one of")
  expect_error(genexpo(TRUE, c(TRUE, FALSE)), "^genexpo: ")
})


test_that("indgap()'s box marks the pooled row and reads the ratios", {
  own <- data.frame(
    group = c("A", "B", "A and B", "Ref"),
    pct = c(31.2, 27.5, 30.8, 18.9),
    lower = c(28, 15, 28.5, 18.4),
    upper = c(34, 40, 33, 19.4),
    stringsAsFactors = FALSE
  )
  g <- indgap(
    own,
    "group",
    "pct",
    "lower",
    "upper",
    ref = "Ref",
    pooled = "A and B"
  )
  txt <- ls_flat(g)
  expect_match(txt, "A and B (pooled)", fixed = TRUE)
  expect_match(
    txt,
    "Above the reference (interval above 1): A and A and B",
    fixed = TRUE
  )
  expect_match(
    txt,
    "Not distinguishable from 1 (interval includes it): B",
    fixed = TRUE
  )
  low <- own
  low$pct <- c(10, 11, 10.5, 18.9)
  low$lower <- low$pct - 1
  low$upper <- low$pct + 1
  g2 <- indgap(low, "group", "pct", "lower", "upper", ref = "Ref")
  txt2 <- ls_flat(g2)
  expect_match(txt2, "Below the reference (interval below 1)", fixed = TRUE)
  expect_false(grepl("Above the reference", txt2, fixed = TRUE))
})


test_that("exptrnd()'s reading follows the direction of the odds", {
  up <- ls_flat(exptrnd(cases = c(60, 95, 70), n = c(500, 520, 300)))
  expect_match(up, "The odds rise with exposure", fixed = TRUE)
  expect_match(up, "At the highest level (2)", fixed = TRUE)
  flat <- ls_flat(exptrnd(cases = c(50, 51, 49), n = c(500, 500, 500)))
  expect_match(flat, "finds no steady change", fixed = TRUE)
  down <- ls_flat(exptrnd(cases = c(120, 60, 20), n = c(500, 500, 500)))
  expect_match(down, "The odds fall with exposure", fixed = TRUE)
})


test_that("exptrnd() refuses a constant outcome", {
  expect_error(exptrnd(c(0, 0, 0), c(10, 10, 10)), "no trend to test")
  expect_error(exptrnd(c(10, 10, 10), c(10, 10, 10)), "no trend to test")
})

test_that("unknown quality flags warn and count as ok", {
  d <- data.frame(
    g = c("A", "B", "C"),
    e = c(10, 12, 11),
    lo = c(8, 10, 9),
    hi = c(12, 14, 13),
    q = c("", "zz", NA),
    stringsAsFactors = FALSE
  )
  expect_warning(
    {
      r <- indgap(d, "g", "e", "lo", "hi", ref = "C", qual = "q")
    },
    "\"zz\" not recognised"
  )
  expect_identical(r$usewithcau, c(FALSE, FALSE))
  d$q <- c("x", "..", "")
  expect_silent(indgap(d, "g", "e", "lo", "hi", ref = "C", qual = "q"))
})

test_that("the results print, with the notes on how to read them", {
  d <- data.frame(
    REF_DATE = "2019/2022",
    GEO = "Canada",
    `Indigenous identity` = c(
      "First Nations",
      "Inuk (Inuit)",
      "Non-Indigenous"
    ),
    VALUE = c("45.8", "", "61.5"),
    STATUS = c("", "x", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  t <- sctidy(d)
  expect_output(expect_invisible(print(t)), "Statistics Canada table")
  expect_match(ls_flat(t), "ok 2, x 1", fixed = TRUE)
  expect_match(ls_flat(t), "x suppressed", fixed = TRUE)
  expect_match(ls_flat(t[c("GEO", "VALUE")]), "GEO VALUE", fixed = TRUE)
  g <- indgap(
    .ih()[
      .ih()$period == "2019/2022" &
        .ih()$age == "Total, 18 years and over" &
        .ih()$indicator == "Diabetes",
    ],
    "group",
    "estimate",
    "lower",
    "upper",
    qual = "quality"
  )
  expect_output(print(g), "Gap against a reference group")
  expect_match(ls_flat(g), "Inuit *", fixed = TRUE)
  expect_match(ls_flat(g), "use with caution", fixed = TRUE)
  expect_match(
    ls_flat(exptrnd(c(60, 95, 70), c(500, 520, 300))),
    "The odds rise with exposure",
    fixed = TRUE
  )
  expect_match(
    ls_flat(exptrnd(c(50, 51, 49), c(500, 500, 500))),
    "finds no steady change",
    fixed = TRUE
  )
  expect_match(
    ls_flat(exptrnd(c(120, 60, 20), c(500, 500, 500))),
    "The odds fall with exposure",
    fixed = TRUE
  )
  expect_match(
    ls_flat(exptrnd(c(0, 5), c(10, 10))),
    "Haldane-Anscombe",
    fixed = TRUE
  )
  expect_match(ls_flat(trcprog()), "13 of 94", fixed = TRUE)
  expect_match(
    ls_flat(trcprog(done = c(0, 0), yearend = c(2019, 2020))),
    "never",
    fixed = TRUE
  )
  expect_match(ls_flat(trcprog()[c("rate", "years")]), "rate", fixed = TRUE)
  l <- limfrom(c(10, 20), se = c(1, 2))
  expect_s3_class(l, c("limfrom", "data.frame"))
  expect_output(print(l), "Intervals from published precision")
  expect_match(ls_flat(l), "[8.04, 11.96]", fixed = TRUE)
  expect_match(ls_flat(l["estimate"]), "estimate", fixed = TRUE)
  low <- data.frame(
    group = c("A", "B", "Ref"),
    pct = c(10, 11, 18.9),
    lower = c(9, 10, 18.4),
    upper = c(11, 12, 19.4),
    stringsAsFactors = FALSE
  )
  expect_match(
    ls_flat(indgap(low, "group", "pct", "lower", "upper", ref = "Ref")),
    "Below the reference (interval below 1)",
    fixed = TRUE
  )
  expect_match(ls_flat(g[c("group", "ratio")]), "ratio", fixed = TRUE)
})

test_that("a reference that looks like a total is said to be one", {
  d <- data.frame(
    g = c("A", "B", "Total"),
    e = c(10, 12, 11),
    lo = c(8, 10, 9),
    hi = c(12, 14, 13),
    stringsAsFactors = FALSE
  )
  expect_warning(
    indgap(d, "g", "e", "lo", "hi", ref = "Total"),
    "looks like a total"
  )
})
