# Alcohol and drug use (R/adx.R): the blood alcohol concentration against
# the formulas written out and the values the sources print, the flags of
# the Canadian surveys, and the scores of the screening questionnaires.

test_that("ebac() is Widmark's formula, in either units", {
  # Matthews and Miller (1979): 5.14 oz / (r lb) - 0.015 t, with 0.6 fl oz
  # of ethanol in a 14 g drink
  mm <- function(d, lb, t, r) max(0, 5.14 * 0.6 * d / (lb * r) - 0.015 * t)
  expect_equal(
    ebac(4, 180, 2, unit = "lb"),
    mm(4, 180, 2, 0.73),
    tolerance = 2e-3
  )
  expect_equal(
    ebac(3, 130, 1, sex = "female", unit = "lb"),
    mm(3, 130, 1, 0.66),
    tolerance = 2e-3
  )
  # grams over (10 r kg), less beta t
  expect_identical(ebac(4, 80, 2), 14 * 4 / (10 * 0.73 * 80) - 0.015 * 2)
  expect_identical(
    ebac(4, 80, 2, r = 0.68, beta = 0.02, grams = 10),
    10 * 4 / (10 * 0.68 * 80) - 0.02 * 2
  )
  # pounds are kilograms times 0.45359237
  expect_identical(ebac(4, 180, 2, unit = "lb"), ebac(4, 180 * 0.45359237, 2))
  # never below 0, vectors recycle from length 1, a missing input is NA
  expect_identical(ebac(1, 80, 10), 0)
  expect_identical(ebac(4, 80, 0:3), vapply(0:3, function(t) ebac(4, 80, t), 0))
  expect_identical(ebac(c(2, NA), 80, 1), c(ebac(2, 80, 1), NA))
  expect_identical(
    ebac(4, 80, 2, sex = c("male", "female")),
    c(ebac(4, 80, 2), ebac(4, 80, 2, r = 0.66))
  )
  expect_identical(ebac(4, 80), ebac(4, 80, 0))
})

test_that("ebac() with a height uses the factors of Seidl, Jensen and Alt", {
  rf <- 0.31223 - 0.006446 * 62 + 0.004466 * 164
  rm <- 0.31608 - 0.004821 * 75 + 0.004632 * 178
  expect_identical(
    ebac(5, 62, 3, sex = "female", height = 164, beta = 0.017, grams = 13.6),
    13.6 * 5 / (10 * rf * 62) - 0.017 * 3
  )
  expect_identical(
    ebac(4, 75, 3, sex = "male", height = 178, beta = 0.0017, grams = 13.6),
    13.6 * 4 / (10 * rm * 75) - 0.0017 * 3
  )
  # the survey's public file: four drinks in three hours give medians of
  # 0.13 (women) and 0.09 (men), at the rate its guide writes
  expect_equal(
    round(
      ebac(4, 62, 3, sex = "female", height = 164, beta = 0.0017, grams = 13.6),
      2
    ),
    0.13
  )
  expect_equal(
    round(
      ebac(4, 75, 3, sex = "male", height = 178, beta = 0.0017, grams = 13.6),
      2
    ),
    0.09
  )
  # one height per person, with the weight in pounds, and r in its place
  expect_identical(
    ebac(4, c(62, 75), 3, sex = c("female", "male"), height = c(164, 178)),
    c(ebac(4, 62, 3, r = rf), ebac(4, 75, 3, r = rm))
  )
  expect_identical(
    ebac(4, 62 / 0.45359237, 3, sex = "female", height = 164, unit = "lb"),
    ebac(4, (62 / 0.45359237) * 0.45359237, 3, sex = "female", height = 164)
  )
  expect_identical(
    ebac(4, 62, 3, height = 164, r = 0.6),
    ebac(4, 62, 3, r = 0.6)
  )
})

test_that("ebac() refuses bad input", {
  errs(
    "ebac",
    list(
      list("`drinks` must be numeric", "a", 80),
      list("`weight` must be numeric", 4, NULL),
      list("`hours` must be numeric", 4, 80, "x"),
      list("`drinks` must be numbers from 0", -1, 80),
      list("`weight` must be positive numbers", 4, 0),
      list("`hours` must be numbers from 0", 4, 80, -1),
      list("`sex` must be \"male\" or \"female\"", 4, 80, 1, sex = "other"),
      list(
        "`sex` must be \"male\" or \"female\"",
        4,
        80,
        1,
        sex = NA_character_
      ),
      list("`height` must be numeric", 4, 80, 1, height = "tall"),
      list("`height` must be positive numbers", 4, 80, 1, height = 0),
      list(
        "the weight and height give a factor that is not positive",
        4,
        400,
        1,
        sex = "female",
        height = 100
      ),
      list(
        "the arguments must have one length, or length 1",
        4,
        c(60, 70),
        1,
        sex = c("male", "female", "male"),
        height = 170
      ),
      list("`r` must be numeric", 4, 80, 1, r = "a"),
      list("`r` must be positive numbers", 4, 80, 1, r = 0),
      list("`beta` must be one number from 0", 4, 80, 1, beta = -0.1),
      list("`beta` must be one number from 0", 4, 80, 1, beta = c(0.1, 0.2)),
      list("`grams` must be one positive number", 4, 80, 1, grams = 0),
      list(
        "the arguments must have one length, or length 1",
        c(1, 2),
        c(60, 70, 80)
      ),
      list("`unit` must be one of", 4, 80, unit = "stone")
    )
  )
})

test_that("stddrk() counts standard drinks by a country's grams", {
  # Canada: 341 mL of 5% beer, 142 mL of 12% wine and 43 mL of 40% spirits
  # are one drink of 13.45 g each
  expect_equal(
    stddrk(c(341, 142, 43), c(5, 12, 40)),
    c(1, 1, 1),
    tolerance = 0.01
  )
  expect_identical(stddrk(750, 13), 750 * 13 / 100 * 0.789 / 13.45)
  expect_identical(stddrk(750, 13, "us"), 750 * 0.13 * 0.789 / 14)
  expect_identical(stddrk(750, 13, "uk"), 750 * 0.13 * 0.789 / 8)
  expect_identical(stddrk(750, 13, "au"), stddrk(750, 13, "who"))
  expect_identical(stddrk(750, 13, grams = 13.6), 750 * 0.13 * 0.789 / 13.6)
  expect_identical(stddrk(c(100, NA), 10), c(stddrk(100, 10), NA))
  expect_identical(stddrk(0, 40), 0)
  errs(
    "stddrk",
    list(
      list("`vol` must be numeric", "a", 5),
      list("`abv` must be numeric", 341, NULL),
      list("`vol` must be numbers from 0", -1, 5),
      list("`abv` must be per cent, from 0 to 100", 341, 500),
      list("`abv` must be per cent, from 0 to 100", 341, -1),
      list("`grams` must be one positive number", 341, 5, grams = 0),
      list("`grams` must be one positive number", 341, 5, grams = "a"),
      list(
        "the arguments must have one length, or length 1",
        c(1, 2),
        c(1, 2, 3)
      ),
      list("`std` must be one of", 341, 5, std = "fr")
    )
  )
})

test_that("hvydrk() and lrdg() are the survey's two flags", {
  d <- c(2, 3, 4, 5, NA)
  # four or more for women, five or more for men
  expect_identical(hvydrk(d, "female"), c(FALSE, FALSE, TRUE, TRUE, NA))
  expect_identical(hvydrk(d, "male"), c(FALSE, FALSE, FALSE, TRUE, NA))
  # more than two for women, more than three for men
  expect_identical(lrdg(d, "female"), c(FALSE, TRUE, TRUE, TRUE, NA))
  expect_identical(lrdg(d, "male"), c(FALSE, FALSE, TRUE, TRUE, NA))
  expect_identical(
    lrdg(d, c("female", "male", "male", "female", "male")),
    c(FALSE, FALSE, TRUE, TRUE, NA)
  )
  expect_identical(hvydrk(6, c("female", "male"), f = 7, m = 6), c(FALSE, TRUE))
  expect_identical(lrdg(3, "male", m = 2), TRUE)
  for (f in c("hvydrk", "lrdg")) {
    errs(
      f,
      list(
        list("`drinks` must be numeric", "a", "male"),
        list("`drinks` must be numbers from 0", -1, "male"),
        list("`sex` must be \"male\" or \"female\"", 1, "x"),
        list("`sex` must be \"male\" or \"female\"", 1, 1),
        list("`f` must be one number from 0", 1, "male", f = -1),
        list("`m` must be one number from 0", 1, "male", m = c(1, 2)),
        list(
          "the arguments must have one length, or length 1",
          c(1, 2, 3),
          c("male", "female")
        )
      )
    )
  }
})

a <- rbind(
  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  c(3, 2, 2, 1, 1, 0, 1, 1, 0, 2),
  c(4, 4, 4, 3, 3, 2, 3, 3, 4, 4)
)

test_that("auditc() scores the AUDIT and the AUDIT-C", {
  r <- auditc(a, sex = c("female", "male", "male"))
  expect_named(r, c("total", "auditc", "zone", "cpos"))
  expect_identical(r[["total"]], c(1, 13, 34))
  expect_identical(r[["auditc"]], c(1, 7, 12))
  expect_identical(r[["zone"]], c("low", "hazardous", "dependence"))
  expect_identical(r[["cpos"]], c(FALSE, TRUE, TRUE))
  # the zones change at 8, 16 and 20
  z <- function(t) {
    m <- matrix(0, 1, 10)
    m[1, seq_len(t %/% 4)] <- 4
    m[1, t %/% 4 + 1] <- t %% 4
    auditc(m)[["zone"]]
  }
  expect_identical(
    vapply(c(7, 8, 15, 16, 19, 20), z, ""),
    c("low", "hazardous", "hazardous", "harmful", "harmful", "dependence")
  )
  # the AUDIT-C alone: positive from 3 for women and from 4 for men
  c3 <- auditc(rbind(c(1, 1, 1), c(2, 1, 1)), sex = "female")
  expect_identical(c3[["total"]], c(NA_real_, NA_real_))
  expect_identical(c3[["zone"]], c(NA_character_, NA_character_))
  expect_identical(c3[["cpos"]], c(TRUE, TRUE))
  expect_identical(
    auditc(rbind(c(1, 1, 1), c(2, 1, 1)), sex = "male")[["cpos"]],
    c(FALSE, TRUE)
  )
  expect_identical(auditc(a)[["cpos"]], c(NA, NA, NA))
  # a data frame, and a missing answer
  expect_identical(auditc(as.data.frame(a)), auditc(a))
  b <- a
  b[2, 5] <- NA
  expect_identical(auditc(b)[["total"]], c(1, NA, 34))
  expect_identical(auditc(b)[["zone"]], c("low", NA, "dependence"))
  errs(
    "auditc",
    list(
      list(
        "`x` must be a numeric matrix or data frame of 10 columns",
        a[, 1:9]
      ),
      list("`x` must be a numeric matrix or data frame", "a"),
      list("`x` must be a numeric matrix or data frame", data.frame(x = "a")),
      list("item 10 must be scored 0, 2 or 4", cbind(a[, 1:9], 3)),
      list("item 2 must be scored 0 to 4", cbind(1, 5, a[, 3:10])),
      list(
        "`sex` must be \"male\" or \"female\", one or one per row",
        a,
        sex = "x"
      ),
      list(
        "`sex` must be \"male\" or \"female\", one or one per row",
        a,
        sex = c("male", "female")
      )
    )
  )
})

test_that("dast10() scores the DAST-10, item 3 reversed", {
  r <- dast10(rbind(
    c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(1, 1, 0, 1, 1, 0, 1, 1, 0, 1),
    c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  ))
  expect_identical(r[["total"]], c(1, 8, 0, 10))
  expect_identical(r[["level"]], c("low", "substantial", "none", "severe"))
  # item 3 answered yes counts nothing; t of the other nine answered yes
  lv <- function(t) {
    m <- matrix(0, 1, 10)
    m[1, 3] <- 1
    k <- c(1, 2, 4:10)
    m[1, k] <- as.numeric(seq_along(k) <= t)
    dast10(m)[["level"]]
  }
  expect_identical(
    vapply(c(0, 1, 2, 3, 5, 6, 8, 9), lv, ""),
    c(
      "none",
      "low",
      "low",
      "moderate",
      "moderate",
      "substantial",
      "substantial",
      "severe"
    )
  )
  expect_identical(dast10(matrix(c(1, NA, rep(0, 8)), 1))[["total"]], NA_real_)
  errs(
    "dast10",
    list(
      list(
        "`x` must be a numeric matrix or data frame of 10 columns",
        matrix(0, 1, 9)
      ),
      list(
        "item 4 must be 0 .no. or 1 .yes.",
        matrix(c(0, 0, 0, 2, rep(0, 6)), 1)
      )
    )
  )
})

test_that("cudit() scores the CUDIT-R", {
  r <- cudit(rbind(
    c(1, 0, 0, 0, 0, 0, 0, 0),
    c(4, 3, 2, 1, 2, 1, 2, 4),
    c(2, 2, 2, 1, 0, 0, 0, 0),
    c(4, 4, 0, 0, 0, 0, 0, 4)
  ))
  expect_identical(r[["total"]], c(1, 19, 7, 12))
  expect_identical(r[["level"]], c("low", "disorder", "low", "disorder"))
  expect_identical(cudit(matrix(c(4, 4, rep(0, 6)), 1))[["level"]], "hazardous")
  expect_identical(cudit(matrix(c(NA, rep(0, 7)), 1))[["level"]], NA_character_)
  errs(
    "cudit",
    list(
      list(
        "`x` must be a numeric matrix or data frame of 8 columns",
        matrix(0, 1, 10)
      ),
      list("item 8 must be scored 0, 2 or 4", matrix(c(rep(0, 7), 1), 1)),
      list("item 1 must be scored 0 to 4", matrix(c(9, rep(0, 7)), 1))
    )
  )
})

test_that("canast() weights the ASSIST items as the surveys derive the score", {
  # the user guide: CAA00 0/2/3/4/6, CAA01 0/3/4/5/6, CAA02 0/4/5/6/7,
  # CAA03 0/5/6/7/8, CAA04 and CAA05 0/3/6
  w <- list(
    c(0, 2, 3, 4, 6),
    c(0, 3, 4, 5, 6),
    c(0, 4, 5, 6, 7),
    c(0, 5, 6, 7, 8),
    c(0, 3, 6),
    c(0, 3, 6)
  )
  for (j in 1:6) {
    for (k in seq_along(w[[j]])) {
      m <- matrix(0, 1, 6)
      m[1, j] <- k - 1
      expect_identical(canast(m)[["total"]], w[[j]][k])
    }
  }
  r <- canast(rbind(
    c(1, 0, 0, 0, 0, 0),
    c(3, 2, 0, 0, 0, 0),
    c(4, 4, 3, 3, 2, 2),
    c(1, 1, 0, 0, 0, 0),
    c(4, 4, 4, 2, 2, 0),
    c(0, 0, 0, 0, 0, 0)
  ))
  expect_identical(r[["total"]], c(2, 8, 37, 5, 31, 0))
  # 0 to 3 low, 4 to 26 moderate, 27 to 39 high
  expect_identical(
    r[["level"]],
    c("low", "moderate", "high", "moderate", "high", "low")
  )
  expect_identical(canast(matrix(c(4, 4, 4, 4, 2, 2), 1))[["total"]], 39)
  expect_identical(canast(matrix(c(1, NA, 0, 0, 0, 0), 1))[["total"]], NA_real_)
  errs(
    "canast",
    list(
      list(
        "`x` must be a numeric matrix or data frame of 6 columns",
        matrix(0, 1, 8)
      ),
      list("item 5 must be answered 0 to 2", matrix(c(1, 0, 0, 0, 3, 0), 1)),
      list("item 2 must be answered 0 to 4", matrix(c(1, 5, 0, 0, 0, 0), 1)),
      list("item 1 must be answered 0 to 4", matrix(c(1.5, 0, 0, 0, 0, 0), 1))
    )
  )
})
