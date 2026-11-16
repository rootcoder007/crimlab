# The census profile and the Calls to Action 19 gaps, on the made-up tables
# of helper-igh.R: the source package's tests of these functions.

test_that(
  paste0(
    "cnsprof() reads every group, the ",
    "hierarchy, the sample and the symbols"
  ),
  {
    raw <- .ipp()
    p <- cnsprof(.ipp())
    expect_s3_class(p, "cnsprof")
    expect_identical(nrow(p), nrow(raw) * 9L)
    expect_setequal(
      unique(p$group),
      c(
        "Total",
        "Indigenous identity",
        "Single Indigenous responses",
        "First Nations",
        "M\u00e9tis",
        "Inuit",
        "Multiple Indigenous responses",
        "Indigenous responses not included elsewhere",
        "Non-Indigenous"
      )
    )
    fc <- unique(p[
      p$characteristic == "Foster children",
      c("level", "parid", "rootid", "sample")
    ])
    expect_identical(c(fc$level, fc$parid, fc$rootid), c(2L, 37L, 13L))
    expect_identical(fc$sample, "25%")
    expect_identical(unique(p$sample[p$chaid == 923L]), "25%")
    # one value against the published cell
    r <- raw[
      raw$GEO_NAME == "Canada" &
        raw$CHARACTERISTIC_ID == "38" &
        raw$AGE_GROUP == "0 to 14 years" &
        raw$RESIDENCE == "On reserve",
    ]
    v <- p$value[
      p$geo == "Canada" &
        p$chaid == 38L &
        p$age == "0 to 14 years" &
        p$residence == "On reserve" &
        p$group == "First Nations"
    ]
    expect_identical(v, as.numeric(r$FIRST_NATIONS))
    # a suppressed or unavailable cell has no value
    sx <- raw$SYMBOL_INUK %in% c("x", "F", "..", "...")
    if (any(sx)) {
      i <- which(sx)[1]
      pv <- p$value[
        p$geo == raw$GEO_NAME[i] &
          p$chaid == as.integer(raw$CHARACTERISTIC_ID[i]) &
          p$age == raw$AGE_GROUP[i] &
          p$residence == raw$RESIDENCE[i] &
          p$group == "Inuit"
      ]
      expect_true(is.na(pv))
    }
    expect_error(
      cnsprof(data.frame(a = 1)),
      "not an Indigenous Population Profile table"
    )
    expect_error(cnsprof(1), "must be a data frame")
    bad <- raw[1:2, c("GEO_NAME", "CHARACTERISTIC_ID", "CHARACTERISTIC_NAME")]
    expect_error(cnsprof(bad), "no group columns")
  }
)


test_that(
  paste0(
    "profshr(): the published counts, and the interval of a 25% sample with ",
    "random rounding"
  ),
  {
    raw <- .ipp()
    p <- cnsprof(.ipp())
    fc <- profshr(p, "Foster children", "Total - Children in census families")
    r <- fc[
      fc$geo == "Canada" &
        fc$age == "0 to 14 years" &
        fc$residence == "Total - Residence on or off reserve" &
        fc$group == "First Nations",
    ]
    sel <- function(id) {
      raw[
        raw$GEO_NAME == "Canada" &
          raw$CHARACTERISTIC_ID == id &
          raw$AGE_GROUP == "0 to 14 years" &
          raw$RESIDENCE == "Total - Residence on or off reserve",
        "FIRST_NATIONS"
      ]
    }
    n <- as.numeric(sel("38"))
    d <- as.numeric(sel("20"))
    expect_identical(c(r$numerator, r$denominator), c(n, d))
    q <- n / d
    expect_identical(r$estimate, 100 * q)
    v <- 0.75 * q * (1 - q) / (0.25 * d) + 4 * (1 + q^2) / d^2
    ne <- q * (1 - q) / v
    z <- stats::qnorm(0.975)
    wl <- (q + z^2 / (2 * ne) - z * sqrt(q * (1 - q) / ne + z^2 / (4 * ne^2))) /
      (1 + z^2 / ne)
    expect_equal(r$lower, 100 * wl, tolerance = testthat::testthat_tolerance())
    # the default denominator is the parent: persons not in census families
    dflt <- profshr(p, 38)
    expect_identical(
      attr(dflt, "characteristics")[["denominator"]],
      "Total - Persons not in census families"
    )
    # a design effect widens it
    w2 <- profshr(
      p,
      "Foster children",
      "Total - Children in census families",
      deff = 2
    )
    expect_gt(w2$upper[1] - w2$lower[1], fc$upper[1] - fc$lower[1])
    expect_error(profshr(p, "No such thing"), "names no characteristic")
    expect_error(profshr(p, "Participation rate"), "has no parent")
    expect_error(profshr(data.frame(), 1), "must come from cnsprof")
    expect_error(profshr(p, c(1, 2)), "one characteristic id or name")
  }
)


test_that("cnscwb(): the four components worked from the published figures", {
  raw <- .ipp()
  p <- cnsprof(.ipp())
  h <- cnsprof(.ihh())
  cw <- cnscwb(p, h)
  fn <- cw[cw$geo == "Canada" & cw$group == "First Nations", ]
  tot <- "Total - Residence on or off reserve"
  g <- function(id, age) {
    as.numeric(raw[
      raw$GEO_NAME == "Canada" &
        raw$CHARACTERISTIC_ID == id &
        raw$AGE_GROUP == age &
        raw$RESIDENCE == tot,
      "FIRST_NATIONS"
    ])
  }
  inc <- g("807", "Total - Age groups") *
    g("809", "Total - Age groups") /
    g("1", "Total - Age groups")
  expect_identical(fn$incpercap, inc)
  expect_identical(
    fn$income,
    100 * (log(inc) - log(2650)) / (log(75000) - log(2650))
  )
  expect_identical(
    fn$highsch25p,
    100 * (1 - g("546", "25 years and over") / g("545", "25 years and over"))
  )
  expect_identical(
    fn$labour,
    (g("629", "25 to 64 years") + 100 - g("631", "25 to 64 years")) / 2
  )
  hr <- .ihh()
  hv <- function(id) {
    as.numeric(hr[
      hr$GEO_NAME == "Canada" & hr$CHARACTERISTIC_ID == id,
      "IND_HHLD (2)"
    ])
  }
  expect_identical(
    fn$housing,
    (100 * hv("38") / hv("37") + 100 * (1 - hv("54") / hv("52"))) / 2
  )
  expect_equal(
    fn$cwb,
    mean(c(fn$income, fn$education, fn$housing, fn$labour)),
    tolerance = 1e-12
  )
  # fewer than 250 people: no score
  small <- cw[is.finite(cw$population) & cw$population < 250, ]
  expect_gt(nrow(small), 0)
  expect_true(all(is.na(small$cwb)))
  # without the household table, no housing component and no score
  nh <- cnscwb(p)
  expect_true(all(is.na(nh$housing)))
  expect_true(all(is.na(nh$cwb)))
  expect_false(all(is.na(nh$income)))
  expect_error(cnscwb(data.frame()), "must come from cnsprof")
  expect_error(
    cnscwb(p, dwell = data.frame()),
    "`dwell` must come from cnsprof"
  )
})


test_that("profshr() names `num` and `den` as the signature does", {
  p <- cnsprof(.ipp())
  expect_error(profshr(p, c(1, 2)), "`num` must be one characteristic")
  expect_error(profshr(p, "No such thing"), "`num` \\(No such thing\\) names")
  expect_error(profshr(p, 38, den = NA), "`den` must be one characteristic")
})


test_that("cnsprof()'s box keeps residence, age, gender and the hierarchy", {
  p <- cnsprof(.ipp())
  fc <- p[p$characteristic == "Foster children" & p$geo == "Canada", ]
  txt <- ls_flat(fc)
  for (k in c(
    "residence",
    "On reserve",
    "Off reserve",
    "0 to 14 years",
    "Total - Gender",
    "parid",
    "rootid",
    "dguid",
    "geolvl",
    "tnrlf"
  )) {
    expect_match(txt, k, fixed = TRUE)
  }
  expect_match(txt, "None carries a quality flag.", fixed = TRUE)
  fc$quality[1L] <- "E"
  expect_match(ls_flat(fc), "1 carries a quality flag (E 1)", fixed = TRUE)
})


test_that("profshr()'s box keeps gender, age, DGUID and quality", {
  p <- cnsprof(.ipp())
  fc <- profshr(p, "Foster children", "Total - Children in census families")
  fc <- fc[fc$age == "0 to 14 years" & fc$group == "First Nations", ]
  txt <- ls_flat(utils::head(fc, 6L))
  for (k in c(
    "gender",
    "Total - Gender",
    "dguid",
    "2021A000011124",
    "quality"
  )) {
    expect_match(txt, k, fixed = TRUE)
  }
  expect_match(txt, "Highest share: Canada, Off reserve", fixed = TRUE)
})

test_that(
  paste0(
    "trcind() covers every Call 19 indicator, each with a direction or a ",
    "reason for none"
  ),
  {
    ti <- trcind()
    c19 <- c(
      "infant mortality",
      "maternal health",
      "suicide",
      "mental health",
      "addictions",
      "life expectancy",
      "birth rates",
      "infant and child health issues",
      "chronic diseases",
      "illness and injury incidence",
      "availability of appropriate health services"
    )
    expect_setequal(unique(ti$indicator[ti$call == "19"]), c19)
    expect_true(all(ti$better %in% c("higher", "lower", "")))
    expect_true(all(nzchar(ti$note[!nzchar(ti$better)])))
    expect_setequal(
      unique(ti$indicator[ti$call == "55"]),
      c("children in care", "educational attainment", "income attainment")
    )
  }
)


test_that("trchlthgap() takes tables of your own and their all-ages label", {
  expect_identical(
    trchlthgap(
      .ih(),
      .ile(),
      .iim(),
      .isu(),
      agetot = c("Total, 18 years and over", "Total, 15 years and over")
    ),
    trchlthgap(.ih(), .ile(), .iim(), .isu())
  )
  le <- data.frame(
    period = "2017/2019",
    group = c("First Nations", "Inuit", "Non-Indigenous"),
    estimate = c(74.2, 70.9, 81.8),
    lower = c(73.6, 69.1, 81.6),
    upper = c(74.8, 72.7, 82.0),
    quality = "ok",
    stringsAsFactors = FALSE
  )
  g <- trchlthgap(health = NULL, lifeexp = le, infmort = NULL, suiide = NULL)
  expect_s3_class(g, "trchlthgap")
  expect_equal(g$disratio, 81.8 / c(74.2, 70.9), tolerance = 1e-12)
  # a health table with its own age label, or none
  h <- .ih()
  h <- h[h$age == "Total, 18 years and over", ]
  h2 <- h
  h2$age <- "All adults"
  a <- trchlthgap(health = h, lifeexp = NULL, infmort = NULL, suiide = NULL)
  b <- trchlthgap(
    health = h2,
    lifeexp = NULL,
    infmort = NULL,
    suiide = NULL,
    agetot = "All adults"
  )
  expect_equal(a, b, tolerance = 1e-12)
  h2$age <- NULL
  expect_equal(
    trchlthgap(health = h2, lifeexp = NULL, infmort = NULL, suiide = NULL),
    a,
    tolerance = 1e-12
  )
  expect_error(trchlthgap(agetot = 1), "trchlthgap: `agetot` must be")
  expect_error(trchlthgap(agetot = letters[1:3]), "trchlthgap: `agetot`")
})


test_that("the suicidal-thoughts gap uses the latest period only", {
  d <- data.frame(
    measure = "m",
    group = c("A", "B", "A", "B"),
    period = c("2015", "2015", "2019", "2019"),
    estimate = c(10, 5, 12, 6),
    lower = c(8, 4, 10, 5),
    upper = c(12, 6, 14, 7),
    quality = "A",
    stringsAsFactors = FALSE
  )
  d$measure <- "Seriously contemplated suicide in the past 12 months"
  r <- trchlthgap(suiide = d, ref = "B")
  expect_true(all(r$period == "2019"))
  expect_identical(r$estimate[r$group == "A"], 12)
})


test_that("cnsstd(): the census age structure, as weights", {
  p <- cnsprof(.ipp())
  ages <- c("0 to 14 years", "25 years and over")
  s <- cnsstd(p, ages = ages)
  expect_equal(sum(s$weight), 1, tolerance = 1e-15)
  expect_identical(nrow(s), 2L)
  raw <- .ipp()
  r <- raw[
    raw$GEO_NAME == "Canada" &
      raw$CHARACTERISTIC_ID == "1" &
      raw$AGE_GROUP == "25 years and over" &
      raw$RESIDENCE == "Total - Residence on or off reserve",
    "IND_IDENT (2)"
  ]
  expect_identical(s$population[s$age == "25 years and over"], as.numeric(r))
  expect_error(cnsstd(p), "no population for: 0 to 4 years")
  expect_error(cnsstd(p, grp = "Nobody"), "no population for group")
  expect_error(cnsstd(p, ages = NA_character_), "`ages` must name age groups")
  # it feeds the standardisation directly
  sp <- stdprev(c(10, 20), c(8, 18), c(12, 22), std = s$population)
  expect_s3_class(sp, "stdprev")
  two <- p
  two$dguid[two$geo == "Lakeside"] <- "other"
  two$geo[two$geo == "Lakeside"] <- "Canada"
  expect_error(cnsstd(two, ages = ages), "names 2 places")
})

test_that("the census results print, with the notes on reading them", {
  p <- cnsprof(.ipp())
  fc <- p[p$characteristic == "Foster children" & p$geo == "Canada", ]
  txt <- ls_flat(fc)
  for (k in c("residence", "On reserve", "0 to 14 years", "Total - Gender")) {
    expect_match(txt, k, fixed = TRUE)
  }
  expect_match(txt, "None carries a quality flag.", fixed = TRUE)
  expect_match(ls_flat(p), "carry a quality flag", fixed = TRUE)
  fc$quality[1L] <- "E"
  expect_match(ls_flat(fc), "1 carries a quality flag", fixed = TRUE)
  expect_output(print(p), "Indigenous Population Profile")
  expect_match(ls_flat(p["geo"]), "geo", fixed = TRUE)
  s <- profshr(p, "Foster children", "Total - Children in census families")
  expect_output(print(s), "Share of a population")
  expect_match(ls_flat(s), "a lower bound on the uncertainty", fixed = TRUE)
  expect_match(ls_flat(s["geo"]), "geo", fixed = TRUE)
  cw <- cnscwb(p, cnsprof(.ihh()))
  expect_output(print(cw), "Community Well-Being index, from the census")
  expect_match(ls_flat(cw), "Never a measure of wellness", fixed = TRUE)
  expect_match(
    ls_flat(cnscwb(p)),
    "Housing needs the household table",
    fixed = TRUE
  )
  expect_match(ls_flat(cw["geo"]), "geo", fixed = TRUE)
  g <- trchlthgap(.ih(), .ile(), .iim(), .isu())
  expect_output(print(g), "Call 19: health gaps")
  expect_match(ls_flat(g), "Life expectancy at birth", fixed = TRUE)
  expect_match(ls_flat(g["group"]), "group", fixed = TRUE)
})

test_that(
  paste0(
    "trchlthgap(): each gap is indgap's ratio, turned so that above 1 is ",
    "worse"
  ),
  {
    g <- trchlthgap(.ih(), .ile(), .iim(), .isu())
    expect_s3_class(g, "trchlthgap")
    h <- .ih()
    d <- h[
      h$period == "2019/2022" &
        h$age == "Total, 18 years and over" &
        h$indicator == "Mood disorder",
    ]
    ig <- indgap(d, "group", "estimate", "lower", "upper", qual = "quality")
    r <- g[g$measure == "Mood disorder" & g$group == "First Nations", ]
    expect_identical(r$disratio, ig$ratio[ig$group == "First Nations"])
    expect_identical(r$lower, ig$ratiolo[ig$group == "First Nations"])
    # higher is better: the ratio is inverted, and so is its interval
    le <- .ile()
    le <- le[le$period == "2016/2021", ]
    il <- indgap(le, "group", "estimate", "lower", "upper", qual = "quality")
    r2 <- g[g$measure == "Life expectancy at birth" & g$group == "Inuit", ]
    expect_identical(r2$disratio, 1 / il$ratio[il$group == "Inuit"])
    expect_identical(r2$lower, 1 / il$ratioup[il$group == "Inuit"])
    expect_identical(r2$period, "2016/2021")
    expect_true(r2$worse)
    # the latest period with values: 2014/2016, not the unavailable 1994/1996
    expect_identical(
      unique(g$period[g$measure == "Infant mortality rate"]),
      "2014/2016"
    )
    expect_identical(g$worse, is.finite(g$lower) & g$lower > 1)
    # the suicidal-thoughts table has no non-Indigenous rows, so suicidal
    # thoughts are not among the gaps
    expect_false(
      "Seriously contemplated suicide in the past 12 months" %in% g$measure
    )
    only <- trchlthgap(lifeexp = .ile())
    expect_identical(unique(only$measure), "Life expectancy at birth")
    expect_error(trchlthgap(NULL, NULL, NULL, NULL), "no table given measures")
    expect_error(trchlthgap(health = 1), "`health` must be a data frame")
    expect_error(trchlthgap(lifeexp = 1), "must be data frames or NULL")
    expect_error(trchlthgap(suiide = 1), "`suiide` must be a data frame")
    su <- .isu()
    su <- rbind(
      su,
      transform(
        su[1, ],
        group = "Non-Indigenous",
        estimate = 2,
        lower = 1.8,
        upper = 2.2
      )
    )
    gs <- trchlthgap(suiide = su)
    expect_identical(unique(gs$period), "2017")
    su$period <- NULL
    expect_identical(unique(trchlthgap(suiide = su)$period), "2017")
  }
)

test_that("a cnsprof object with its columns missing is refused, not read", {
  # the class is trusted until a column is read without a further check;
  # the other arguments are checked first (see test-igr.R)
  cp <- structure(data.frame(), class = c("cnsprof", "data.frame"))
  expect_error(cnsstd(cp), "must come from cnsprof()", fixed = TRUE)
})
