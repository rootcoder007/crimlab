# Mental health and wellbeing from published tables: the source package's
# tests of these functions.

test_that(
  paste0(
    "scests() puts each estimate beside its limits, cell by cell and ",
    "measure by measure"
  ),
  {
    d <- data.frame(
      REF_DATE = "2022",
      GEO = "Canada",
      `Indigenous identity` = rep(c("First Nations ", "Metis"), each = 5),
      Statistics = rep(
        c(
          "Percent",
          "Low 95% confidence interval",
          "High 95% confidence interval",
          "Number of persons",
          "Low 95% confidence interval, number of persons"
        ),
        2
      ),
      UOM = rep(c("Percent", "Percent", "Percent", "Persons", "Persons"), 2),
      VALUE = c(47.1, 45.2, 49.0, 1000, 900, 50.3, 48.1, 52.5, 800, 700),
      STATUS = c("", "", "", "", "", "F", "", "", "", ""),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    e <- scests(d)
    p <- e[e$measure == "Percent", ]
    expect_identical(p[["Indigenous identity"]], c("First Nations", "Metis"))
    expect_identical(p$estimate, c(47.1, NA))
    expect_identical(p$lower, c(45.2, NA))
    expect_identical(p$upper, c(49.0, NA))
    expect_identical(p$cisou, c("published", NA))
    expect_identical(as.character(p$quality), c("ok", "F"))
    n <- e[e$measure == "Number of persons", ]
    expect_identical(n$estimate, c(1000, 800))
    # a lower limit published without its upper is kept, and the pair is not
    # called published
    expect_identical(n$lower, c(900, 700))
    expect_identical(n$upper, c(NA_real_, NA_real_))
    expect_identical(n$cisou, c(NA_character_, NA_character_))
    expect_identical(attr(e, "conflvl"), 0.95)
    d3 <- d[c(1:3, 1), ]
    expect_error(scests(d3), "two estimates for one cell")
    # one place name, two geographies: told apart by DGUID, not merged
    dg <- data.frame(
      GEO = "Division No.  1",
      DGUID = c("2021A00031001", "2021A00034601"),
      Statistics = "Count",
      VALUE = c(10, 20),
      STATUS = "",
      stringsAsFactors = FALSE
    )
    eg <- scests(dg)
    expect_identical(eg$DGUID, c("2021A00031001", "2021A00034601"))
    expect_identical(eg$estimate, c(10, 20))
  }
)


test_that(
  paste0(
    "scests(): coefficients of variation, standard errors, and limits named ",
    "in every form"
  ),
  {
    z <- stats::qnorm(0.975)
    cv <- scests(data.frame(
      GEO = "Canada",
      Characteristics = c("Rate", "Coefficient of variation for rate"),
      VALUE = c(12.5, 8),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_equal(cv$lower, 12.5 * exp(-z * sqrt(log(1 + 0.08^2))),
      tolerance = 1e-12)
    expect_equal(cv$upper, 12.5 * exp(z * sqrt(log(1 + 0.08^2))),
      tolerance = 1e-12)
    expect_identical(cv$cisou, "cv")
    se <- scests(data.frame(
      GEO = "A",
      Statistics = c("Average score", "Standard error"),
      VALUE = c(250, 2),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(se$upper, 250 + z * 2)
    expect_identical(se$cisou, "se")
    ct <- scests(data.frame(
      GEO = "A",
      Statistics = c(
        "Count",
        "Confidence interval lower bound of count",
        "Confidence interval upper bound of count"
      ),
      VALUE = c(100, 80, 120),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(c(ct$lower, ct$upper), c(80, 120))
    yr <- scests(data.frame(
      GEO = "A",
      Statistics = c(
        "2015",
        "2020",
        "95% confidence interval lower bound, 2015",
        "95% confidence interval upper bound, 2020"
      ),
      VALUE = c(10, 20, 8, 23),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(yr$lower, c(8, NA))
    expect_identical(yr$upper, c(NA, 23))
    lim <- scests(data.frame(
      GEO = "A",
      Statistics = c(
        "Percentage",
        "95% confidence interval (lower limit)",
        "95% confidence interval (upper limit)"
      ),
      VALUE = c(30, 25, 35),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(c(lim$lower, lim$upper), c(25, 35))
    # a name that fits three measures: the one whose estimates the limits
    # contain
    py <- data.frame(
      GEO = rep(c("A", "B"), each = 5),
      Characteristics = rep(
        c(
          "Number of potential years of life lost",
          "Crude rate of potential years of life lost",
          "Age-standardized rate of potential years of life lost",
          "Low 95% confidence interval, potential years of life lost",
          "High 95% confidence interval, potential years of life lost"
        ),
        2
      ),
      VALUE = c(5000, 900, 1100, 1000, 1200, 7000, 800, 950, 870, 1030),
      STATUS = "",
      stringsAsFactors = FALSE
    )
    pe <- scests(py)
    asr <- pe[
      pe$measure == "Age-standardized rate of potential years of life lost",
    ]
    expect_identical(asr$lower, c(1000, 870))
    expect_true(all(is.na(pe$lower[
      pe$measure != "Age-standardized rate of potential years of life lost"
    ])))
    expect_identical(attr(pe, "unmatched"), character())
    # one cell in 200 whose published limits miss its estimate: the match holds,
    # and the cell is marked
    many <- do.call(
      rbind,
      lapply(seq_len(200), function(i) {
        transform(py[1:5, ], GEO = paste0("R", i))
      })
    )
    many$VALUE[many$GEO == "R7"][3] <- 1300
    me <- scests(many)
    ma <- me[
      me$measure == "Age-standardized rate of potential years of life lost",
    ]
    expect_identical(attr(me, "unmatched"), character())
    expect_identical(ma$cisou[ma$GEO == "R7"], "inconsistent")
    expect_identical(sum(ma$cisou == "published"), 199L)
    # a limit no single measure fits is said, not guessed
    amb <- scests(data.frame(
      GEO = "A",
      Statistics = c("Men", "Women", "Low 95% confidence interval"),
      VALUE = c(10, 12, 5),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(attr(amb, "unmatched"), "Low 95% confidence interval")
    # no statistics column: every value its own estimate, with no limits
    nc <- scests(data.frame(
      GEO = c("A", "B"),
      Gender = "Total",
      UOM = "Number",
      VALUE = c(5, 6),
      STATUS = "",
      stringsAsFactors = FALSE
    ))
    expect_identical(nc$estimate, c(5, 6))
    expect_identical(nc$measure, c("Number", "Number"))
    expect_true(all(is.na(nc$lower)))
  }
)

# shares and Wald limits that imply an effective sample size of exactly n

test_that(
  paste0(
    "ordcmp(): the probability of a higher rating, brute force over every ",
    "pair"
  ),
  {
    lv <- c("poor", "fair", "good", "excellent")
    p1 <- c(0.1, 0.2, 0.3, 0.4)
    p0 <- c(0.25, 0.25, 0.3, 0.2)
    w1 <- .wald(p1, 400)
    w0 <- .wald(p0, 900)
    d <- data.frame(
      g = rep(c("A", "R"), each = 4),
      cat = rep(lv, 2),
      pct = 100 * c(p1, p0),
      lo = c(w1$lo, w0$lo),
      hi = c(w1$hi, w0$hi),
      stringsAsFactors = FALSE
    )
    d <- rbind(
      d,
      data.frame(
        g = "A",
        cat = "not specified",
        pct = 3,
        lo = 2,
        hi = 4,
        stringsAsFactors = FALSE
      )
    )
    r <- ordcmp(d, "g", "cat", "pct", "lo", "hi", lvls = lv, ref = "R")
    x <- rep(1:4, 100 * p1)
    y <- rep(1:4, 100 * p0)
    cmp <- outer(x, y, function(a, b) (a > b) + 0.5 * (a == b))
    expect_equal(
      r$probbet,
      mean(cmp),
      tolerance = testthat::testthat_tolerance()
    )
    expect_equal(
      r$delta,
      2 * mean(cmp) - 1,
      tolerance = testthat::testthat_tolerance()
    )
    expect_equal(r$neff, 400, tolerance = testthat::testthat_tolerance())
    expect_equal(r$neffref, 900, tolerance = testthat::testthat_tolerance())
    # the U-statistic variance from the placements of each sample in the other
    v <- mean((rowMeans(cmp) - mean(cmp))^2) /
      400 +
      mean((colMeans(cmp) - mean(cmp))^2) / 900
    sl <- sqrt(v) / (mean(cmp) * (1 - mean(cmp)))
    expect_equal(
      r$upper,
      stats::plogis(stats::qlogis(mean(cmp)) + stats::qnorm(0.975) * sl),
      tolerance = testthat::testthat_tolerance()
    )
    expect_identical(dim(attr(r, "shares")), c(2L, 4L))
    # identical distributions: one half, interval about it
    d2 <- d[d$cat %in% lv, ]
    d2[d2$g == "A", c("pct", "lo", "hi")] <- d2[
      d2$g == "R",
      c("pct", "lo", "hi")
    ]
    r2 <- ordcmp(d2, "g", "cat", "pct", "lo", "hi", lvls = lv, ref = "R")
    expect_equal(r2$probbet, 0.5)
    expect_lt(r2$lower, 0.5)
    expect_gt(r2$upper, 0.5)
  }
)


test_that("ordcmp(): suppressed categories, caution and bad arguments", {
  lv <- c("low", "high")
  d <- data.frame(
    g = rep(c("A", "B", "R"), each = 2),
    cat = rep(lv, 3),
    pct = c(40, 60, 30, 70, 50, 50),
    lo = c(35, 55, 25, 65, 45, 45),
    hi = c(45, 65, 35, 75, 55, 55),
    q = c("", "", "F", "", "", "E"),
    stringsAsFactors = FALSE
  )
  r <- ordcmp(
    d,
    "g",
    "cat",
    "pct",
    "lo",
    "hi",
    lvls = lv,
    ref = "R",
    qual = "q"
  )
  expect_true(is.na(r$probbet[r$group == "B"]))
  expect_true(all(r$usewithcau))
  expect_error(
    ordcmp(d, "g", "cat", "pct", "lo", "hi", lvls = "low", ref = "R"),
    "two or more"
  )
  expect_error(
    ordcmp(d, "g", "cat", "pct", "lo", "hi", lvls = lv, ref = "Z"),
    "is not a group"
  )
  expect_error(
    ordcmp(d, "g", "cat", "pct", "lo", "hi", lvls = c("x", "y"), ref = "R"),
    "no row of `data`"
  )
  d$q[5] <- "x"
  expect_error(
    ordcmp(d, "g", "cat", "pct", "lo", "hi", lvls = lv, ref = "R", qual = "q"),
    "reference group's distribution"
  )
  expect_error(
    ordcmp(
      rbind(d, d[1, ]),
      "g",
      "cat",
      "pct",
      "lo",
      "hi",
      lvls = lv,
      ref = "A"
    ),
    "has a category twice"
  )
})


test_that("stdprev(): weighted mean and the root of the weighted variances", {
  z <- stats::qnorm(0.975)
  r <- stdprev(c(6, 4, 2), c(5, 3, 1), c(7, 5, 3), std = c(2, 1, 1))
  expect_identical(r$estimate, sum(c(0.5, 0.25, 0.25) * c(6, 4, 2)))
  expect_identical(r$se, sqrt(sum(c(0.5, 0.25, 0.25)^2 * (1 / z)^2)))
  expect_identical(r$upper, r$estimate + z * r$se)
  r2 <- stdprev(
    c(6, 4, 2, 5, 5, 5),
    c(5, 3, 1, 4, 4, 4),
    c(7, 5, 3, 6, 6, 6),
    std = c(1, 1, 1),
    grp = rep(c("a", "b"), each = 3),
    qual = c("", "x", "", "E", "", "")
  )
  expect_true(is.na(r2$estimate[1]))
  expect_true(r2$usewithcau[2])
  expect_identical(r2$estimate[2], 5)
  expect_warning(
    {
      r3 <- stdprev(1, 2, 3, 1)
    },
    "do not contain their estimate"
  )
  expect_true(is.na(r3$estimate))
  expect_error(
    stdprev(c(1, 2), c(0, 1), c(2, 3), c(1, 1, 1)),
    "has 2 age groups"
  )
  expect_error(stdprev(c(1, 2), c(0, 1), 3, c(1, 1)), "same length")
  expect_error(stdprev(1, 0, 2, 1, qual = c("", "")), "one flag per row")
})


test_that(
  paste0(
    "welprof(): wellness first, totals left ",
    "out, suppressed values dropped"
  ),
  {
    w <- welprof(
      .wb(),
      "group",
      "indicator",
      "pct",
      "lo",
      "hi",
      well = c("Satisfied", "Hopeful"),
      distress = "Distress",
      domain = c(
        Hopeful = "spiritual",
        Satisfied = "mental",
        Distress = "mental"
      ),
      qual = "q"
    )
    expect_s3_class(w, "welprof")
    expect_identical(
      w$indicator,
      rep(c("Satisfied", "Hopeful", "Distress"), each = 2)
    )
    expect_identical(
      w$side,
      rep(c("wellness", "wellness", "distress"), each = 2)
    )
    expect_identical(w$domain[3], "spiritual")
    expect_true(w$usewithcau[1])
    expect_error(
      welprof(
        .wb(),
        "group",
        "indicator",
        "pct",
        "lo",
        "hi",
        well = "Nope",
        distress = "Distress"
      ),
      "not an indicator"
    )
    expect_error(
      welprof(
        .wb(),
        "group",
        "indicator",
        "pct",
        "lo",
        "hi",
        well = "Hopeful",
        distress = "Hopeful"
      ),
      "both one of wellness"
    )
    expect_error(
      welprof(
        .wb(),
        "group",
        "indicator",
        "pct",
        "lo",
        "hi",
        well = "Hopeful",
        distress = "Distress",
        domain = "x"
      ),
      "named character vector"
    )
  }
)


test_that("dualcont(): the strength ratio and its log-normal interval", {
  w <- welprof(
    .wb(),
    "group",
    "indicator",
    "pct",
    "lo",
    "hi",
    well = "Hopeful",
    distress = "Distress"
  )
  dc <- dualcont(w, "Hopeful", "Distress")
  z <- stats::qnorm(0.975)
  sl <- sqrt(
    ((log(65) - log(55)) / (2 * z))^2 + ((log(23) - log(17)) / (2 * z))^2
  )
  expect_identical(dc$strratio[1], 3)
  expect_identical(dc$ratiolo[1], exp(log(3) - z * sl))
  expect_error(dualcont(w, "Distress", "Distress"), "not a wellness indicator")
  expect_error(dualcont(w, "Hopeful", "Hopeful"), "not a distress indicator")
  expect_error(dualcont(data.frame(), "a", "b"), "must come from welprof")
  # a profile whose last column was renamed: it used to be read anyway,
  # and the missing column crashed R
  v <- w
  names(v)[length(v)] <- names(v)[1]
  expect_error(dualcont(v, "Hopeful", "Distress"),
               "must come from welprof()", fixed = TRUE)
})


test_that(
  paste0(
    "holbal(): one when every domain is ",
    "level; Atkinson's form otherwise"
  ),
  {
    d <- data.frame(
      g = rep(c("A", "B"), each = 4),
      dom = rep(c("physical", "emotional", "mental", "spiritual"), 2),
      pct = c(50, 50, 50, 50, 80, 80, 20, 20),
      lo = c(47, 47, 47, 47, 77, 77, 17, 17),
      hi = c(53, 53, 53, 53, 83, 83, 23, 23),
      stringsAsFactors = FALSE
    )
    r <- holbal(d, "g", "dom", "pct", "lo", "hi", reps = 2000)
    expect_identical(r$balance[1], 1)
    expect_equal(r$balance[2], sqrt(0.8 * 0.2) / 0.5, tolerance = 1e-12)
    expect_identical(r$level, c(50, 50))
    expect_lte(r$balup[1], 1)
    expect_lt(r$balup[2], 1)
    expect_identical(r$weadom[2], "mental")
    expect_identical(
      holbal(d, "g", "dom", "pct", "lo", "hi", reps = 500, seed = 4)$ballo,
      holbal(d, "g", "dom", "pct", "lo", "hi", reps = 500, seed = 4)$ballo
    )
    wt <- c(physical = 1, emotional = 1, mental = 1, spiritual = 3)
    r2 <- holbal(d, "g", "dom", "pct", "lo", "hi", w = wt, reps = 500)
    w <- wt / 6
    expect_identical(
      r2$balance[2],
      exp(sum(w * log(c(0.8, 0.8, 0.2, 0.2)))) / sum(w * c(0.8, 0.8, 0.2, 0.2))
    )
    expect_identical(r2$level[2], 100 * sum(w * c(0.8, 0.8, 0.2, 0.2)))
    d$q <- c("", "", "", "x", "", "", "", "")
    expect_true(is.na(holbal(
      d,
      "g",
      "dom",
      "pct",
      "lo",
      "hi",
      qual = "q",
      reps = 200
    )$balance[1]))
    expect_error(
      holbal(d, "g", "dom", "pct", "lo", "hi", w = c(a = 1)),
      "named by the domains"
    )
    expect_error(
      holbal(d[d$dom == "mental", ], "g", "dom", "pct", "lo", "hi"),
      "two domains"
    )
    d$pct[1] <- 100
    d$hi[1] <- 100
    expect_error(
      holbal(d, "g", "dom", "pct", "lo", "hi"),
      "strictly between 0 and 100"
    )
  }
)


test_that("mental-health errors name the user's arguments", {
  expect_error(
    welprof(
      .wb(),
      "group",
      "indicator",
      "pct",
      "lo",
      "hi",
      well = mean,
      distress = "Distress"
    ),
    "`well` must be data"
  )
  expect_error(stdprev(1, 0, 2, 1, grp = mean), "`grp` must be data")
  expect_error(stdprev(1, 0, 2, 1, qual = mean), "`qual` must be data")
  d <- data.frame(
    g = rep(c("A", "B"), each = 4),
    dom = rep(c("physical", "emotional", "mental", "spiritual"), 2),
    pct = c(50, 50, 50, 50, 80, 80, 20, 20),
    lo = c(47, 47, 47, 47, 77, 77, 17, 17),
    hi = c(53, 53, 53, 53, 83, 83, 23, 23),
    stringsAsFactors = FALSE
  )
  expect_error(
    holbal(d, "g", "dom", "pct", "lo", "hi", w = mean),
    "`w` must be data"
  )
  oc <- data.frame(
    g = rep(c("A", "R"), each = 2),
    cat = rep(c("low", "high"), 2),
    pct = c(40, 60, 50, 50),
    lo = c(35, 55, 45, 45),
    hi = c(45, 65, 55, 55),
    stringsAsFactors = FALSE
  )
  expect_error(
    ordcmp(oc, "g", "cat", "pct", "lo", "hi", lvls = mean, ref = "R"),
    "`lvls` must be data"
  )
  expect_error(
    ordcmp(
      oc,
      "g",
      "cat",
      "pct",
      "lo",
      "hi",
      lvls = c("low", "high"),
      ref = "R",
      pooled = "Totl"
    ),
    "^ordcmp: `pooled` names \"Totl\""
  )
})


test_that("the help pages' own tables run through the estimators", {
  own <- data.frame(
    group = rep(c("Nation A", "Everyone else"), each = 4),
    rating = rep(c("Poor", "Fair", "Good", "Excellent"), 2),
    pct = c(10, 25, 40, 25, 6, 20, 44, 30),
    lo = c(7, 21, 35, 21, 5, 18, 42, 28),
    hi = c(13, 29, 45, 29, 7, 22, 46, 32),
    stringsAsFactors = FALSE
  )
  oc <- ordcmp(
    own,
    "group",
    "rating",
    "pct",
    "lo",
    "hi",
    lvls = c("Poor", "Fair", "Good", "Excellent"),
    ref = "Everyone else"
  )
  expect_s3_class(oc, "ordcmp")
  sp <- stdprev(
    est = c(12, 9, 6, 8, 6, 4),
    lower = c(9, 6, 3, 7, 5, 3),
    upper = c(15, 12, 9, 9, 7, 5),
    std = c(500, 300, 200),
    grp = rep(c("Nation A", "Everyone else"), each = 3)
  )
  expect_s3_class(sp, "stdprev")
  expect_equal(sp$estimate[sp$group == "Nation A"], 9.9, tolerance = 1e-12)
  wp <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 4),
    item = rep(
      c("Hopeful", "Satisfied", "Rarely hopeful", "Low satisfaction"),
      2
    ),
    pct = c(64, 71, 12, 18, 58, 66, 16, 22),
    lo = c(59, 66, 9, 14, 52, 60, 12, 17),
    hi = c(69, 76, 15, 22, 64, 72, 20, 27),
    stringsAsFactors = FALSE
  )
  w <- welprof(
    wp,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    well = c("Hopeful", "Satisfied"),
    distress = c("Rarely hopeful", "Low satisfaction")
  )
  expect_s3_class(w, "welprof")
  expect_s3_class(frmaud(w), "frmaud")
  expect_s3_class(dualcont(w, "Hopeful", "Rarely hopeful"), "dualcont")
  hb <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 4),
    domain = rep(c("body", "mind", "spirit", "family"), 2),
    pct = c(62, 58, 70, 75, 55, 61, 66, 80),
    lo = c(57, 53, 65, 70, 49, 55, 60, 75),
    hi = c(67, 63, 75, 80, 61, 67, 72, 85),
    stringsAsFactors = FALSE
  )
  expect_s3_class(
    holbal(hb, "community", "domain", "pct", "lo", "hi", reps = 2000),
    "holbal"
  )
})


test_that("stdprev()'s reading compares the groups' prevalences", {
  sp <- stdprev(
    c(5, 6, 1, 2),
    c(4, 5, 0.5, 1.5),
    c(6, 7, 1.5, 2.5),
    std = c(1, 1),
    grp = c("x", "x", "y", "y")
  )
  txt <- ls_flat(sp)
  expect_match(txt, "Highest prevalence: x", fixed = TRUE)
  expect_match(txt, "lowest: y", fixed = TRUE)
  expect_match(txt, "do not overlap", fixed = TRUE)
  sp2 <- stdprev(
    c(5, 6, 5.5, 6.2),
    c(3, 4, 3.5, 4.2),
    c(7, 8, 7.5, 8.2),
    std = c(1, 1),
    grp = c("x", "x", "y", "y")
  )
  expect_match(ls_flat(sp2), "Highest prevalence: y", fixed = TRUE)
  expect_match(ls_flat(sp2), "Their intervals overlap", fixed = TRUE)
})


test_that("the mental health results print, distinctions kept", {
  own <- data.frame(
    group = rep(c("Nation A", "Everyone else", "Total"), each = 4),
    rating = rep(c("Poor", "Fair", "Good", "Excellent"), 3),
    pct = c(10, 25, 40, 25, 6, 20, 44, 30, 8, 22, 42, 28),
    lo = c(7, 21, 35, 21, 5, 18, 42, 28, 6, 20, 40, 26),
    hi = c(13, 29, 45, 29, 7, 22, 46, 32, 10, 24, 44, 30),
    stringsAsFactors = FALSE
  )
  oc <- ordcmp(
    own,
    "group",
    "rating",
    "pct",
    "lo",
    "hi",
    lvls = c("Poor", "Fair", "Good", "Excellent"),
    ref = "Everyone else"
  )
  expect_output(print(oc), "Ordinal comparison")
  expect_match(ls_flat(oc), "Total (pooled)", fixed = TRUE)
  expect_match(ls_flat(oc), "Rate lower than Everyone else", fixed = TRUE)
  expect_match(ls_flat(oc["group"]), "group", fixed = TRUE)
  sp <- stdprev(
    c(12, 9, 6, 8, 6, 4),
    c(9, 6, 3, 7, 5, 3),
    c(15, 12, 9, 9, 7, 5),
    c(500, 300, 200),
    grp = rep(c("A", "B"), each = 3),
    qual = c("E", "", "", "", "", "")
  )
  expect_output(print(sp), "Age-standardised prevalence")
  expect_match(ls_flat(sp), "Highest prevalence: A *", fixed = TRUE)
  expect_match(ls_flat(sp["group"]), "group", fixed = TRUE)
  own4 <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 4),
    item = rep(
      c("Hopeful", "Satisfied", "Rarely hopeful", "Low satisfaction"),
      2
    ),
    pct = c(64, 71, 12, 18, 38, 66, 16, 22),
    lo = c(59, 66, 9, 14, 32, 60, 12, 17),
    hi = c(69, 76, 15, 22, 44, 72, 20, 27),
    stringsAsFactors = FALSE
  )
  w <- welprof(
    own4,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    well = c("Hopeful", "Satisfied"),
    distress = c("Rarely hopeful", "Low satisfaction"),
    domain = c(Hopeful = "spirit")
  )
  expect_output(print(w), "Wellness profile")
  expect_match(ls_flat(w), "Hopeful (wellness)", fixed = TRUE)
  expect_match(ls_flat(w), "spirit", fixed = TRUE)
  expect_match(ls_flat(w), "intervals part on 1 of 4", fixed = TRUE)
  w2 <- welprof(
    own4[own4$community == "Lakeside", ],
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    well = "Hopeful",
    distress = "Rarely hopeful"
  )
  expect_match(ls_flat(w2), "On none of the 2 indicators", fixed = TRUE)
  expect_match(ls_flat(w["group"]), "group", fixed = TRUE)
  dc <- dualcont(w, "Hopeful", "Rarely hopeful")
  expect_output(print(dc), "Dual continuum")
  expect_match(
    ls_flat(dc),
    "2 of 2 groups have the whole interval above 1",
    fixed = TRUE
  )
  expect_match(ls_flat(dc["group"]), "group", fixed = TRUE)
  hb <- holbal(
    data.frame(
      c = rep(c("L", "R"), each = 2),
      d = rep(c("body", "mind"), 2),
      p = c(62, 58, 55, 80),
      lo = c(57, 53, 49, 75),
      hi = c(67, 63, 61, 85),
      stringsAsFactors = FALSE
    ),
    "c",
    "d",
    "p",
    "lo",
    "hi",
    reps = 200
  )
  expect_output(print(hb), "Holistic balance")
  expect_match(ls_flat(hb), "equal, over 2 domains", fixed = TRUE)
  expect_match(ls_flat(hb), "Domain levels (percent)", fixed = TRUE)
  expect_match(ls_flat(hb), "The weakest domain is", fixed = TRUE)
  hw <- holbal(
    data.frame(
      c = "L",
      d = c("body", "mind"),
      p = c(62, 58),
      lo = c(57, 53),
      hi = c(67, 63),
      stringsAsFactors = FALSE
    ),
    "c",
    "d",
    "p",
    "lo",
    "hi",
    reps = 200,
    w = c(body = 2, mind = 1)
  )
  expect_match(ls_flat(hw), "body 0.667", fixed = TRUE)
  attr(hw, "domains") <- NULL
  expect_false(grepl("Domain levels", ls_flat(hw), fixed = TRUE))
  expect_match(ls_flat(hb["group"]), "group", fixed = TRUE)
  s <- scests(data.frame(
    GEO = "Canada",
    Characteristics = c("Rate", "Count"),
    VALUE = c(1, 2),
    stringsAsFactors = FALSE
  ))
  expect_output(print(s), "Statistics Canada table")
})

test_that("the warnings of the table functions start with their name", {
  expect_warning(
    stdprev(c(5, 20), c(8, 18), c(12, 22), std = c(1, 1)),
    "^stdprev: 1 interval\\(s\\) do not contain their estimate"
  )
})
