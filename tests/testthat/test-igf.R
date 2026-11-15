# The frameworks for wellness put to work on published tables: the source
# package's tests of these functions.

test_that(
  paste0(
    "cwbscr(): the four components and their mean, as the methodology sets ",
    "them"
  ),
  {
    r <- cwbscr(20000, 60, 9, 80, 70, 60, 90)
    expect_identical(
      r$income,
      100 * (log(20000) - log(2650)) / (log(75000) - log(2650))
    )
    expect_identical(r$education, 2 / 3 * 60 + 1 / 3 * 9)
    expect_identical(r$housing, 75)
    expect_identical(r$labour, 75)
    expect_identical(r$cwb, mean(c(r$income, r$education, r$housing, r$labour)))
    expect_identical(
      cwbscr(c(1000, 1e6), 50, 5, 50, 50, 50, 50)$income,
      c(0, 100)
    )
    expect_identical(
      cwbscr(20000, 60, 9, 80, 70, 60, 90, bounds = c(2000, 40000))$income,
      100 * (log(20000) - log(2000)) / (log(40000) - log(2000))
    )
    expect_identical(nrow(cwbscr(c(1e4, 2e4, 3e4), 50, 5, 50, 50, 50, 50)), 3L)
    expect_error(cwbscr(1e4, 150, 5, 50, 50, 50, 50), "from 0 to 100")
    expect_error(cwbscr(0, 50, 5, 50, 50, 50, 50), "must be positive")
    expect_error(
      cwbscr(1e4, c(1, 2), c(1, 2, 3), 50, 50, 50, 50),
      "one value per community"
    )
    expect_error(
      cwbscr(1e4, 50, 5, 50, 50, 50, 50, bounds = c(5, 1)),
      "two increasing"
    )
  }
)


test_that(
  paste0(
    "healchg(): the change between periods ",
    "and its interval, wellness first"
  ),
  {
    z <- stats::qnorm(0.975)
    r <- healchg(
      .hc(),
      "g",
      "ind",
      "per",
      "est",
      "lo",
      "hi",
      from = "t1",
      to = "t2",
      well = "well",
      qual = "q"
    )
    expect_identical(r$indicator, c("well", "well", "ill", "ill"))
    expect_identical(r$side, c("wellness", "wellness", "distress", "distress"))
    expect_identical(r$change[1], 5)
    expect_identical(r$upper[1], 5 + z * sqrt(2 * (8 / (2 * z))^2))
    expect_true(is.na(r$change[2]))
    expect_true(r$usewithcau[4])
    expect_true(all(is.na(
      healchg(.hc(), "g", "ind", "per", "est", "lo", "hi", "t1", "t2")$side
    )))
    expect_error(
      healchg(.hc(), "g", "ind", "per", "est", "lo", "hi", "t1", "t9"),
      "is not a period"
    )
    expect_error(
      healchg(
        rbind(.hc(), .hc()),
        "g",
        "ind",
        "per",
        "est",
        "lo",
        "hi",
        "t1",
        "t2"
      ),
      "appear twice"
    )
    d <- .hc()
    d$ind[d$per == "t2"] <- paste0("other", seq_len(sum(d$per == "t2")))
    expect_error(
      healchg(d, "g", "ind", "per", "est", "lo", "hi", "t1", "t2"),
      "in both periods"
    )
  }
)


test_that(
  paste0(
    "commidx(): the weighted sum over the sum of absolute weights, never ",
    "filled"
  ),
  {
    z <- stats::qnorm(0.975)
    d <- data.frame(
      g = rep(c("A", "B"), each = 3),
      ind = rep(c("belong", "stress", "drink"), 2),
      e = c(60, 20, 30, 70, 25, NA),
      lo = c(56, 18, 27, 66, 22, NA),
      hi = c(64, 22, 33, 74, 28, NA),
      q = c("", "", "", "", "", ".."),
      stringsAsFactors = FALSE
    )
    w <- c(belong = 2, stress = -1, drink = -1)
    r <- commidx(d, "g", "ind", "e", "lo", "hi", w = w, qual = "q")
    expect_identical(r$index[1], (2 * 60 - 20 - 30) / 4)
    se <- c(4, 2, 3) / z
    expect_identical(
      r$upper[1],
      (2 * 60 - 20 - 30) / 4 + z * sqrt(sum(w^2 * se^2)) / 4
    )
    expect_true(is.na(r$index[2]))
    expect_identical(attr(r, "weights"), w)
    expect_error(
      commidx(d, "g", "ind", "e", "lo", "hi"),
      "set by the community"
    )
    expect_error(
      commidx(d, "g", "ind", "e", "lo", "hi", w = c(zz = 1)),
      "not an indicator"
    )
    expect_error(
      commidx(rbind(d, d[1, ]), "g", "ind", "e", "lo", "hi", w = w),
      "indicator twice"
    )
  }
)


test_that("mhcont(): proxies per direction, shared and missing ones said", {
  z <- stats::qnorm(0.975)
  d <- data.frame(
    g = "A",
    ind = c("hopeful", "belong", "mp", "family"),
    e = c(50, 60, 40, 70),
    lo = c(46, 56, 36, 66),
    hi = c(54, 64, 44, 74),
    stringsAsFactors = FALSE
  )
  r <- mhcont(
    d,
    "g",
    "ind",
    "e",
    "lo",
    "hi",
    hope = "hopeful",
    belong = c("belong", "family"),
    meaning = "mp",
    purp = "mp"
  )
  expect_s3_class(r, "welprof")
  expect_identical(r$indicator, c("hope", "belonging", "meaning", "purpose"))
  expect_identical(r$estimate[2], 65)
  expect_identical(r$upper[2], 65 + z * sqrt(2 * (8 / (2 * z))^2) / 2)
  expect_identical(r$shared, c(FALSE, FALSE, TRUE, TRUE))
  r2 <- mhcont(d, "g", "ind", "e", "lo", "hi", hope = "hopeful")
  expect_identical(r2$domain[2:4], rep("no published proxy", 3))
  expect_true(all(is.na(r2$estimate[2:4])))
  expect_error(mhcont(d, "g", "ind", "e", "lo", "hi"), "one direction at least")
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = "zz"),
    "not an indicator"
  )
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = 1),
    "must name indicators"
  )
})


test_that("frmaud(): the share of the report that is deficit-framed", {
  d <- data.frame(
    g = rep(c("A", "B"), each = 3),
    ind = rep(c("w", "d1", "d2"), 2),
    e = 1:6,
    lo = 0:5,
    hi = 2:7,
    stringsAsFactors = FALSE
  )
  w <- welprof(
    d[-1, ],
    "g",
    "ind",
    "e",
    "lo",
    "hi",
    well = "w",
    distress = c("d1", "d2")
  )
  f <- frmaud(w)
  expect_identical(f$shrdef, 4 / 5)
  expect_identical(f$defonly, "A")
  expect_output(print(f), "Reported on deficits alone: A")
  expect_output(print(f), "pair each with a strength-framed level")
  expect_error(frmaud(data.frame()), "must come from welprof")
})


test_that("framework errors name the user's arguments", {
  d <- data.frame(
    g = rep(c("A", "B"), each = 2),
    ind = rep(c("hope1", "bel1"), 2),
    e = c(50, 60, 55, 65),
    lo = c(45, 55, 50, 60),
    hi = c(55, 65, 60, 70),
    stringsAsFactors = FALSE
  )
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = "hope1", belong = mean),
    "`belong` must be data"
  )
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = "hope1", purp = mean),
    "`purp` must be data"
  )
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = "hope1", belong = 1),
    "^mhcont: `belong` must name indicators"
  )
  expect_error(
    mhcont(d, "g", "ind", "e", "lo", "hi", hope = "hope1", purp = NA),
    "^mhcont: `purp` must name indicators"
  )
  expect_error(
    commidx(d, "g", "ind", "e", "lo", "hi", w = mean),
    "`w` must be data"
  )
  expect_error(
    healchg(.hc(), "g", "ind", "per", "est", "lo", "hi", "t1", "t2", mean),
    "`well` must be data"
  )
  # the column checks name the arguments, not long names
  expect_error(
    commidx(d, NULL, "ind", "e", "lo", "hi"),
    "`grp` must name a column"
  )
})


test_that("frmaud() refuses a profile with no indicator", {
  wb <- data.frame(
    group = c("A", "A"),
    indicator = c("Hopeful", "Distress"),
    pct = c(60, 20),
    lo = c(55, 17),
    hi = c(65, 23),
    stringsAsFactors = FALSE
  )
  w <- welprof(
    wb,
    "group",
    "indicator",
    "pct",
    "lo",
    "hi",
    well = "Hopeful",
    distress = "Distress"
  )
  expect_identical(frmaud(w)$shrdef, 0.5)
  empty <- w[0, ]
  class(empty) <- class(w)
  expect_error(frmaud(empty), "^frmaud: `x` has no wellness or distress")
})


test_that("the help pages' own tables run through the framework estimators", {
  own <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 2),
    item = rep(c("Hopeful", "Life has meaning"), 2),
    pct = c(64, 70, 58, 75),
    lo = c(59, 65, 52, 70),
    hi = c(69, 75, 64, 80),
    stringsAsFactors = FALSE
  )
  mc <- mhcont(
    own,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    hope = "Hopeful",
    meaning = "Life has meaning",
    purp = "Life has meaning"
  )
  expect_s3_class(mc, "mhcont")
  pl <- data.frame(
    region = c("North", "South", "East", "West", "Coast", "Interior"),
    language = c(41, 18, 27, 33, 22, 38),
    wellbeing = c(68, 55, 60, 66, 57, 63),
    stringsAsFactors = FALSE
  )
  ca <- connassoc(pl, "region", "language", "wellbeing", reps = 2000)
  expect_s3_class(ca, "connassoc")
  hc <- data.frame(
    community = "Lakeside",
    item = rep(c("Good mental health", "Anxiety"), 2),
    survey = rep(c("2018", "2023"), each = 2),
    pct = c(55, 14, 61, 11),
    lo = c(50, 11, 56, 8),
    hi = c(60, 17, 66, 14),
    stringsAsFactors = FALSE
  )
  h <- healchg(
    hc,
    "community",
    "item",
    "survey",
    "pct",
    "lo",
    "hi",
    from = "2018",
    to = "2023",
    well = "Good mental health",
    distress = "Anxiety"
  )
  expect_s3_class(h, "healchg")
  ci <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 3),
    item = rep(c("Speaks the language", "Has a doctor", "Food insecure"), 2),
    pct = c(38, 72, 21, 25, 64, 30),
    lo = c(33, 67, 17, 20, 58, 25),
    hi = c(43, 77, 25, 30, 70, 35),
    stringsAsFactors = FALSE
  )
  x <- commidx(
    ci,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    w = c("Speaks the language" = 2, "Has a doctor" = 1, "Food insecure" = -1)
  )
  expect_s3_class(x, "commidx")
})


test_that("commidx()'s reading compares the groups' indices", {
  h <- .ih()
  h <- h[h$period == "2019/2022" & h$age == "Total, 18 years and over", ]
  w <- c(2, -1)
  names(w) <- c("Has a regular healthcare provider", "Heavy drinking")
  ci <- commidx(
    h,
    "group",
    "indicator",
    "estimate",
    "lower",
    "upper",
    qual = "quality",
    w = w
  )
  txt <- ls_flat(ci)
  expect_match(txt, "Highest index:", fixed = TRUE)
  expect_match(txt, "lowest:", fixed = TRUE)
  expect_match(txt, "A higher index: more of the protective", fixed = TRUE)
})


test_that(
  paste0(
    "connassoc(): Spearman's rho, its exact p as cor.test gives it, and the ",
    "error interval"
  ),
  {
    d <- data.frame(
      p = letters[1:6],
      x = c(10, 14, 12, 19, 25, 22),
      y = c(31, 30, 36, 40, 39, 47),
      xl = c(8, 12, 10, 17, 23, 20),
      xh = c(12, 16, 14, 21, 27, 24),
      yl = c(29, 28, 34, 38, 37, 45),
      yh = c(33, 32, 38, 42, 41, 49)
    )
    r <- connassoc(d, "p", "x", "y", "xl", "xh", "yl", "yh", reps = 2000)
    ct <- suppressWarnings(stats::cor.test(
      d$x,
      d$y,
      method = "spearman",
      exact = TRUE
    ))
    expect_equal(r$rho, unname(ct$estimate), tolerance = 1e-14)
    expect_equal(r$p_value, ct$p.value, tolerance = 1e-8)
    expect_true(r$exact)
    expect_lt(r$rholo, r$rho)
    expect_gt(r$rhoup, r$rho)
    expect_output(print(r), "does not show that one causes the other")
    r2 <- connassoc(d, "p", "x", "y")
    expect_true(is.na(r2$rholo))
    expect_output(print(r2), "exact, every permutation")
    big <- data.frame(
      p = 1:12,
      x = 1:12,
      y = c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11)
    )
    r3 <- connassoc(big, "p", "x", "y", reps = 2000)
    expect_false(r3$exact)
    expect_lt(r3$p_value, 0.01)
    expect_output(print(r3), "permutation)")
    expect_error(connassoc(d[1:2, ], "p", "x", "y"), "three places")
    expect_error(connassoc(d, "p", "x", "y", "xl"), "all four limit columns")
    flat <- data.frame(p = 1:4, x = c(1, 2, 3, 4), y = c(5, 5, 5, 5))
    expect_match(
      ls_flat(connassoc(flat, "p", "x", "y")),
      "undefined",
      fixed = TRUE
    )
    weak <- data.frame(p = 1:5, x = 1:5, y = c(3, 1, 4, 5, 2))
    expect_match(
      ls_flat(connassoc(weak, "p", "x", "y")),
      "weakly",
      fixed = TRUE
    )
    opp <- data.frame(p = 1:5, x = 1:5, y = c(9, 7, 8, 3, 1))
    expect_match(
      ls_flat(connassoc(opp, "p", "x", "y")),
      "opposite orders",
      fixed = TRUE
    )
    mid <- data.frame(p = 1:6, x = 1:6, y = c(2, 1, 6, 3, 4, 5))
    expect_match(
      ls_flat(connassoc(mid, "p", "x", "y")),
      "moderately",
      fixed = TRUE
    )
  }
)

test_that("the framework results print, distinctions and governance kept", {
  own <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 2),
    item = rep(c("Hopeful", "Life has meaning"), 2),
    pct = c(64, 70, 58, 75),
    lo = c(59, 65, 52, 70),
    hi = c(69, 75, 64, 80),
    q = c("", "E", "", "")
  )
  m <- mhcont(
    own,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    hope = "Hopeful",
    meaning = "Life has meaning",
    purp = "Life has meaning",
    qual = "q"
  )
  expect_output(print(m), "the four directions")
  txt <- ls_flat(m)
  expect_match(txt, "yes (shared)", fixed = TRUE)
  expect_match(txt, "none published", fixed = TRUE)
  expect_match(txt, "No published proxy: belonging", fixed = TRUE)
  expect_match(txt, "use with caution", fixed = TRUE)
  m2 <- mhcont(
    own,
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    hope = "Hopeful",
    belong = "Hopeful",
    meaning = "Life has meaning",
    purp = "Life has meaning"
  )
  expect_match(
    ls_flat(m2),
    "Every direction has a published proxy",
    fixed = TRUE
  )
  expect_match(ls_flat(m[c("group", "estimate")]), "estimate", fixed = TRUE)
  cw <- cwbscr(
    c(18000, 31000),
    c(58, 81),
    c(6, 19),
    c(71, 95),
    c(62, 90),
    c(55, 68),
    c(78, 92)
  )
  expect_output(print(cw), "Community Well-Being scores")
  expect_match(ls_flat(cw[c("cwb")]), "cwb", fixed = TRUE)
  hc <- healchg(
    .hc(),
    "g",
    "ind",
    "per",
    "est",
    "lo",
    "hi",
    from = "t1",
    to = "t2",
    well = "well",
    qual = "q"
  )
  expect_output(print(hc), "Healing as change")
  expect_match(ls_flat(hc), "well (wellness)", fixed = TRUE)
  expect_match(ls_flat(hc), "Wellness:", fixed = TRUE)
  expect_match(ls_flat(hc[c("group", "change")]), "change", fixed = TRUE)
  ci <- commidx(
    .hc()[.hc()$per == "t1", ],
    "g",
    "ind",
    "est",
    "lo",
    "hi",
    w = c(well = 2, ill = -1)
  )
  expect_output(print(ci), "Community-weighted index")
  expect_match(ls_flat(ci), "Highest index:", fixed = TRUE)
  expect_match(ls_flat(ci), "never a ranking", fixed = TRUE)
  expect_match(ls_flat(ci["group"]), "group", fixed = TRUE)
  own4 <- data.frame(
    community = rep(c("Lakeside", "Riverbend"), each = 2),
    item = rep(c("Rarely hopeful", "Low satisfaction"), 2),
    pct = c(12, 18, 16, 22),
    lo = c(9, 14, 12, 17),
    hi = c(15, 22, 20, 27)
  )
  w <- welprof(
    rbind(
      own4,
      data.frame(
        community = "Lakeside",
        item = "Hopeful",
        pct = 60,
        lo = 55,
        hi = 65
      )
    ),
    "community",
    "item",
    "pct",
    "lo",
    "hi",
    well = "Hopeful",
    distress = c("Rarely hopeful", "Low satisfaction")
  )
  fa <- frmaud(w)
  expect_output(print(fa), "Framing audit")
  expect_match(
    ls_flat(fa),
    "Reported on deficits alone: Riverbend",
    fixed = TRUE
  )
  expect_match(
    ls_flat(fa),
    "pair each with a strength-framed level",
    fixed = TRUE
  )
})
