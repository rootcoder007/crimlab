# Review cases of the Indigenous health functions, carried from the source
# package.

test_that("qltyrcd() maps another agency's flags onto ok, E and x", {
  expect_identical(
    qltyrcd(c("", "*", "**", "np", NA), suppr = c("**", "np"), caution = "*"),
    c("ok", "E", "x", "x", "ok")
  )
  expect_identical(qltyrcd(c("E", "F", "..")), c("E", "x", "x"))
  expect_error(
    qltyrcd("a", suppr = "a", caution = "a"),
    "both suppressed and caution"
  )
  expect_error(qltyrcd(list(1)), "vector of flags")
  expect_error(qltyrcd("a", suppr = 1), "must be flags")
})


test_that(
  paste0(
    "limfrom() builds limits from a standard ",
    "error, a CV or a margin of error"
  ),
  {
    z <- stats::qnorm(0.975)
    r <- limfrom(c(40, 10), se = c(2, 1))
    expect_identical(r$lower, c(40, 10) - z * c(2, 1))
    expect_identical(r$upper, c(40, 10) + z * c(2, 1))
    expect_identical(limfrom(40, cv = 5)$se, 2)
    expect_identical(limfrom(40, moe = stats::qnorm(0.95) * 2)$se, 2)
    expect_identical(limfrom(40, moe = z * 2, moelev = 0.95)$se, 2)
    expect_error(limfrom(40), "give one of")
    expect_error(limfrom(40, se = 1, cv = 1), "give one of")
    expect_error(limfrom(c(1, 2, 3), se = c(1, 2)), "one per estimate")
    expect_error(limfrom(40, se = -1), "must not be negative")
  }
)


test_that(
  paste0(
    "indgap() refuses duplicate groups and ",
    "warns on a total as the reference"
  ),
  {
    d <- data.frame(
      g = c("A", "A", "R"),
      est = c(10, 11, 12),
      lo = c(9, 10, 11),
      hi = c(11, 12, 13),
      stringsAsFactors = FALSE
    )
    expect_error(indgap(d, "g", "est", "lo", "hi", ref = "R"), "appears twice")
    d2 <- data.frame(
      g = c("A", "Total, all persons"),
      est = c(10, 12),
      lo = c(9, 11),
      hi = c(11, 13),
      stringsAsFactors = FALSE
    )
    expect_warning(
      indgap(d2, "g", "est", "lo", "hi", ref = "Total, all persons"),
      "total"
    )
    d3 <- data.frame(
      g = c(" A", "B", "R"),
      est = c(10, 11, 12),
      lo = c(9, 10, 11),
      hi = c(11, 12, 13),
      stringsAsFactors = FALSE
    )
    r <- indgap(d3, "g", "est", "lo", "hi", ref = "R", pooled = "B")
    expect_identical(r$group, c("A", "B"))
    expect_identical(r$poogrp, c(FALSE, TRUE))
  }
)


test_that(
  paste0(
    "holbal(): the interval holds the ",
    "estimate, and proportions are refused"
  ),
  {
    d <- data.frame(
      g = rep(c("A", "B"), each = 4),
      dom = rep(c("physical", "emotional", "mental", "spiritual"), 2),
      pct = c(50, 52, 48, 51, 80, 70, 20, 30),
      lo = c(45, 47, 43, 46, 75, 65, 15, 25),
      hi = c(55, 57, 53, 56, 85, 75, 25, 35),
      stringsAsFactors = FALSE
    )
    r <- holbal(d, "g", "dom", "pct", "lo", "hi", reps = 1000, seed = 2)
    expect_true(all(r$ballo <= r$balance & r$balance <= r$balup))
    expect_true(all(r$ballo >= 0 & r$balup <= 1))
    p <- d
    p[c("pct", "lo", "hi")] <- p[c("pct", "lo", "hi")] / 100
    expect_error(holbal(p, "g", "dom", "pct", "lo", "hi"), "percent|proportion")
  }
)


test_that("connassoc(): the Monte Carlo p-value is (b + 1) / (R + 1)", {
  big <- data.frame(p = 1:12, x = 1:12, y = 1:12)
  r <- connassoc(big, "p", "x", "y", reps = 999)
  expect_false(r$exact)
  expect_gt(r$p_value, 0)
  expect_equal(r$p_value * 1000, round(r$p_value * 1000), tolerance = 1e-9)
})


test_that(
  paste0(
    "healchg(): with `distress` given, an ",
    "indicator in neither list is other"
  ),
  {
    d <- data.frame(
      g = "A",
      ind = rep(c("well", "ill", "income"), each = 2),
      per = rep(c("t1", "t2"), 3),
      est = c(50, 55, 20, 18, 30, 31),
      lo = c(46, 51, 18, 16, 28, 29),
      hi = c(54, 59, 22, 20, 32, 33),
      stringsAsFactors = FALSE
    )
    r <- healchg(
      d,
      "g",
      "ind",
      "per",
      "est",
      "lo",
      "hi",
      "t1",
      "t2",
      well = "well",
      distress = "ill"
    )
    expect_identical(
      r$side[match(c("well", "ill", "income"), r$indicator)],
      c("wellness", "distress", "other")
    )
    r2 <- healchg(
      d,
      "g",
      "ind",
      "per",
      "est",
      "lo",
      "hi",
      "t1",
      "t2",
      well = "well"
    )
    expect_identical(r2$side[r2$indicator == "income"], "distress")
    d$lo[1] <- 52
    expect_warning(
      healchg(d, "g", "ind", "per", "est", "lo", "hi", "t1", "t2"),
      "t1"
    )
  }
)


test_that("commidx() keeps its confidence level", {
  d <- data.frame(
    g = rep(c("A", "B"), each = 2),
    ind = rep(c("u", "v"), 2),
    est = c(40, 20, 50, 10),
    lo = c(36, 17, 46, 8),
    hi = c(44, 23, 54, 12),
    stringsAsFactors = FALSE
  )
  ci <- commidx(
    d,
    "g",
    "ind",
    "est",
    "lo",
    "hi",
    w = c(u = 1, v = -1),
    conflvl = 0.9
  )
  expect_identical(attr(ci, "conflvl"), 0.9)
})


test_that("stdprev() caps the limit at `scale`, which may be Inf", {
  s <- stdprev(c(98, 99), c(95, 96), c(100, 100), c(0.5, 0.5))
  expect_lte(s$upper, 100)
  r <- stdprev(c(98, 99), c(95, 96), c(101, 102), c(0.5, 0.5), scale = Inf)
  expect_gt(r$upper, 100)
})


test_that("profshr(): rounding base and the on-reserve full canvass", {
  p <- cnsprof(.ipp())
  a <- profshr(p, "Foster children", "Total - Children in census families")
  b <- profshr(
    p,
    "Foster children",
    "Total - Children in census families",
    roubase = 1
  )
  ok <- is.finite(a$lower) & is.finite(b$lower) & a$estimate > 0
  expect_true(all(
    b$upper[ok] - b$lower[ok] <= a$upper[ok] - a$lower[ok] + 1e-12
  ))
  on <- a$residence == "On reserve" & ok
  if (any(on)) {
    c4 <- profshr(
      p,
      "Foster children",
      "Total - Children in census families",
      smpfrac = 0.5
    )
    expect_identical(c4$lower[on], a$lower[on])
  }
  expect_error(profshr(p, "Foster children", roubase = 0), "roubase")
})


test_that("cnscwb(): minpop and the DGUID column", {
  p <- cnsprof(.ipp())
  h <- cnsprof(.ihh())
  cw <- cnscwb(p, h, minpop = 65)
  expect_true("dguid" %in% names(cw))
  expect_true(all(is.na(cw$cwb[is.finite(cw$population) & cw$population < 65])))
  expect_error(cnscwb(p, h, minpop = -1), "minpop")
})


test_that("trchlthgap(): a different reference, and duplicates warn", {
  hi <- .ih()
  g <- trchlthgap(hi, NULL, NULL, NULL)
  expect_true(all(c("disratio", "lower", "upper") %in% names(g)))
  expect_error(trchlthgap(hi, NULL, NULL, NULL, ref = "Nobody"), "\"Nobody\"")
  dup <- rbind(hi, hi[hi$indicator == "Diabetes", ])
  expect_warning(
    {
      g2 <- trchlthgap(dup, NULL, NULL, NULL)
    },
    "Diabetes .*more than one row per group"
  )
  expect_false("Diabetes" %in% g2$measure)
})


# The source package's coverage tests of these functions.

.cov_tbl <- function() {
  data.frame(
    g = c("A", "A", "B", "B"),
    ind = rep(c("w", "d"), 2),
    e = c(60, 20, 50, 30),
    lo = c(55, 15, 45, 25),
    hi = c(65, 25, 55, 35),
    q = "",
    stringsAsFactors = FALSE
  )
}

# from # from test-coverage-1_0_0.R
test_that(
  paste0(
    "ordcmp() refuses what is not a table, a column it cannot find, and a ",
    "bad quality column"
  ),
  {
    expect_error(
      ordcmp(1, "g", "r", "e", "lo", "hi", lvls = "a", ref = "A"),
      "must be a data frame"
    )
    d <- data.frame(
      g = "A",
      r = "x",
      e = 1,
      lo = 0,
      hi = 2,
      stringsAsFactors = FALSE
    )
    expect_error(
      ordcmp(d, "g", "r", "e", "lo", "nope", lvls = "x", ref = "A"),
      "column `nope` not found"
    )
    expect_error(
      ordcmp(d, "g", "r", "e", "lo", "hi", lvls = "x", ref = "A", qual = "zz"),
      "column `zz` not found"
    )
  }
)

test_that(
  paste0(
    "ordcmp(): an interval that misses its estimate is not compared, and a ",
    "zero-width one gives no interval"
  ),
  {
    lv <- c("Low", "High")
    d <- data.frame(
      g = rep(c("R", "A", "B"), each = 2),
      r = rep(lv, 3),
      e = c(40, 60, 30, 70, 50, 50),
      lo = c(38, 58, 45, 68, 50, 50), # A's first limit lies above its estimate
      # B's limits equal its estimates: no sampling error
      hi = c(42, 62, 32, 72, 50, 50),
      stringsAsFactors = FALSE
    )
    expect_warning(
      {
        oc <- ordcmp(d, "g", "r", "e", "lo", "hi", lvls = lv, ref = "R")
      },
      "does not contain its estimate"
    )
    expect_true(is.na(oc$probbet[oc$group == "A"]))
    b <- oc[oc$group == "B", ]
    expect_false(is.na(b$probbet))
    expect_true(is.na(b$lower))
    expect_true(is.na(b$upper))
  }
)

test_that(
  paste0(
    "stdprev() refuses a bad standard and ",
    "group labels of the wrong length"
  ),
  {
    expect_error(
      stdprev(c(10, 20), c(8, 18), c(12, 22), std = c(0, 0)),
      "non-negative counts"
    )
    expect_error(
      stdprev(c(10, 20), c(8, 18), c(12, 22), std = c(1, NA)),
      "non-negative counts"
    )
    expect_error(
      stdprev(c(10, 20), c(8, 18), c(12, 22), std = c(1, 1), grp = "only one"),
      "one label per row"
    )
  }
)

test_that(
  paste0(
    "welprof() and the column helper refuse a non-table, a missing column ",
    "and non-text indicators"
  ),
  {
    d <- .cov_tbl()
    expect_error(
      welprof(1, "g", "ind", "e", "lo", "hi", well = "w", distress = "d"),
      "must be a data frame"
    )
    expect_error(
      welprof(d, "g", "ind", "e", "lo", "nope", well = "w", distress = "d"),
      "column `nope` not found"
    )
    expect_error(
      welprof(
        d,
        "g",
        "ind",
        "e",
        "lo",
        "hi",
        well = "w",
        distress = "d",
        qual = "zz"
      ),
      "column `zz` not found"
    )
    expect_error(
      welprof(d, "g", "ind", "e", "lo", "hi", well = 1, distress = "d"),
      "must name indicators"
    )
    expect_error(
      welprof(
        d,
        "g",
        "ind",
        "e",
        "lo",
        "hi",
        well = "w",
        distress = NA_character_
      ),
      "must name indicators"
    )
  }
)

test_that(
  paste0(
    "holbal() and connassoc() refuse a non-table and name a missing quality ",
    "column"
  ),
  {
    d <- .cov_tbl()
    expect_error(holbal(1, "g", "ind", "e", "lo", "hi"), "must be a data frame")
    expect_error(
      holbal(d, "g", "ind", "e", "lo", "hi", qual = "zz"),
      "column `zz` not found"
    )
    expect_error(connassoc(1, "p", "x", "y"), "must be a data frame")
  }
)

test_that(
  paste0(
    "scests() picks the Statistics column beside a Characteristics column, ",
    "and skips an empty measure"
  ),
  {
    grp <- rep(c("First Nations", "Inuit"), each = 4)
    stat <- rep(
      c(
        "Rate (A)",
        "Rate (B)",
        "Low 95% confidence interval",
        "High 95% confidence interval"
      ),
      2
    )
    x <- data.frame(
      REF_DATE = "2022",
      GEO = "Canada",
      Characteristics = grp,
      Statistics = stat,
      VALUE = c(50, NA, 45, 55, 40, NA, 35, 45),
      STATUS = "",
      stringsAsFactors = FALSE
    )
    e <- scests(x)
    expect_s3_class(e, "data.frame")
    # the measure with no values cannot be what the limits belong to
    expect_true(all(
      e$lower <= e$estimate & e$estimate <= e$upper,
      na.rm = TRUE
    ))
    expect_true(any(
      e$lower ==
        c(45, 35)[match(e$Characteristics, c("First Nations", "Inuit"))],
      na.rm = TRUE
    ))
  }
)

test_that("trcprog() refuses years out of order", {
  expect_error(trcprog(c(1, 2), c(2017, 2016)), "must be increasing")
})

test_that(
  paste0(
    "cnsstd() refuses a non-profile, non-text ages, and a geography that ",
    "names several places"
  ),
  {
    expect_error(cnsstd(data.frame(a = 1)), "must come from cnsprof")
    p <- cnsprof(.ipp())
    expect_error(cnsstd(p, ages = 1), "must name age groups")
    d <- p[p$geo == "Canada", ]
    d$geo <- "Twice"
    d2 <- d
    d2$dguid <- paste0(d2$dguid, "b")
    pp <- rbind(d, d2)
    class(pp) <- class(p)
    expect_error(cnsstd(pp, geo = "Twice"), "names 2 places")
  }
)

test_that(
  paste0(
    "trchlthgap() adds the suicidal-thoughts gap when the table carries ",
    "reference rows"
  ),
  {
    # reaches trchlthgap's suicidal-thoughts gap, added when the table has the
    # reference group
    s <- .isu()
    r <- s[s$group == "Inuit", ]
    r$group <- "Non-Indigenous"
    r$estimate <- r$estimate / 2
    r$lower <- r$lower / 2
    r$upper <- r$upper / 2
    g <- trchlthgap(
      health = NULL,
      lifeexp = NULL,
      infmort = NULL,
      suiide = rbind(s, r)
    )
    expect_identical(unique(g$indicator), "suicide")
    expect_identical(unique(g$period), "2017")
    inu <- g[g$group == "Inuit", ]
    expect_equal(inu$disratio, 2, tolerance = 1e-8)
  }
)

test_that(
  paste0(
    "re-check of 1.0.0: the battery's last ",
    "raw errors are refused in words"
  ),
  {
    cp <- structure(
      data.frame(),
      class = c("cnsprof", "data.frame")
    )
    expect_error(
      cnscwb(cp, groups = NULL),
      "`groups` must name one group or more"
    )
    expect_error(
      cnscwb(cp, groups = character(0)),
      "`groups` must name one group or more"
    )
    expect_error(
      cnscwb(cp, groups = c("a", "")),
      "`groups` must name one group or more"
    )
  }
)

test_that("the completion-year interval widens with the horizon", {
  w <- vapply(
    c(94, 150, 300),
    function(tot) {
      p <- trcprog(total = tot)
      p$cmplbylate - p$cmplbyearly
    },
    0
  )
  expect_true(all(diff(w) > 0))
})

test_that(
  paste0(
    "sctidy() leaves the labels of a table ",
    "with one dimension as they are"
  ),
  {
    # reaches .sctidy_members returning the table unchanged when fewer
    # than two dimensions
    x <- data.frame(
      REF_DATE = "2021",
      GEO = "Canada",
      Sex = c("M", "M"),
      COORDINATE = c("1.1", "1.1"),
      VALUE = c(1, 2),
      stringsAsFactors = FALSE
    )
    # the geography alone between GEO and the coordinate: one dimension, nothing
    # to tell apart
    one <- data.frame(
      REF_DATE = "2021",
      GEO = "Canada",
      COORDINATE = "1",
      VALUE = 1,
      stringsAsFactors = FALSE
    )
    expect_identical(sctidy(one)$VALUE, 1)
    s <- sctidy(x)
    expect_identical(s$Sex, c("M", "M"))
    expect_identical(s$VALUE, c(1, 2))
  }
)


test_that("the results print in every shape they come in", {
  et <- exptrnd(c(60, 95, 70), c(500, 520, 300))
  expect_output(expect_invisible(print(et)), "Odds by level of exposure")
  expect_match(ls_flat(et["dose"]), "dose", fixed = TRUE)
  p <- trcprog()
  expect_output(expect_invisible(print(p)), "Calls to Action: progress")
  attr(p, "total") <- NULL
  expect_match(ls_flat(p), "13 of 94", fixed = TRUE)
  expect_match(ls_flat(limfrom(1, se = 1)[0, ]), "(no rows)", fixed = TRUE)
  sp <- stdprev(c(10, 10), c(8, 8), c(12, 12), std = 1, grp = c("A", "B"))
  expect_match(ls_flat(sp), "Every group has the same prevalence",
               fixed = TRUE)
})

test_that("arguments that are no data are refused in words", {
  expect_error(limfrom(new.env(), se = 1), "not an environment")
  expect_error(limfrom(1i, se = 1), "not this kind of value")
  expect_error(limfrom(factor("a"), se = 1), "not a factor")
  expect_error(limfrom(list(1), se = 1), "not a list")
  expect_error(limfrom("abc", se = 1), "not text such as \"abc\"")
  expect_error(limfrom(TRUE, se = 1), "not TRUE")
  expect_true(is.na(limfrom(NA, se = 1)$lower))
  expect_error(limfrom(1, se = 1, conflvl = 2), "strictly between 0 and 1")
  expect_error(indgap(data.frame(g = "A", e = 1, l = 0, u = 2), "g", "e",
                      "l", "u", ref = 1), "single non-empty string")
  expect_error(genexpo(TRUE, TRUE, partial = "zz"), "must be one of")
  expect_error(mhcont(data.frame(g = "A"), "g", "zz", "e", "l", "u",
                      hope = "h"), "column `zz` not found")
  long <- strrep("x", 600)
  expect_error(indgap(data.frame(g = "A"), long, "e", "l", "u"), long,
               fixed = TRUE)
  expect_error(polex(1.5), "`seed` must be a single whole number")
  expect_error(profshr(cnsprof(.ipp()), 2.5), "names no characteristic")
})

test_that("a NULL seed draws one, as the session's stream would", {
  set.seed(4)
  a <- polex(NULL)
  set.seed(4)
  b <- polex(NULL)
  expect_identical(a, b)
  d <- data.frame(p = 1:4, x = c(1, 3, 2, 4), y = c(2, 1, 4, 3))
  expect_true(is.finite(connassoc(d, "p", "x", "y", seed = NULL)$rho))
})

test_that("places without both indicators are left out, their rows named", {
  d <- data.frame(p = 1:5, x = c(1, NA, 3, 4, 5), y = c(2, 1, 4, 3, 5))
  r <- connassoc(d, "p", "x", "y")
  expect_identical(rownames(r$data), c("1", "3", "4", "5"))
  expect_identical(r$n, 4L)
})

test_that("tables of odd layouts are read as published", {
  w <- data.frame(GEO = "Canada", `Sex (2):Men[1]` = c(1, 2),
                  check.names = FALSE)
  t <- sctidy(w)
  expect_identical(as.character(t$quality), c("ok", "ok"))
  expect_true(is.na(sctidy(data.frame(GEO = "C", VALUE = "infinity"))$VALUE))
  s <- scests(data.frame(GEO = "C", `Characteristics (3)` = c("A", "B"),
                         VALUE = c(1, 2), check.names = FALSE))
  expect_identical(nrow(s), 2L)
  m <- scests(data.frame(GEO = "C",
                         Measure = c("Rate", "Standard error of rate"),
                         VALUE = c(10, 1)))
  expect_identical(m$cisou, "se")
  l <- scests(data.frame(GEO = "C", L = I(list(1, NULL)),
                         Statistics = c("Rate", "Count"), VALUE = c(1, 2)))
  expect_identical(nrow(l), 2L)
})

test_that("trchlthgap() reads numeric periods and binds unlike ones", {
  le <- .ile()
  le$period <- rep(c(2016, 2021), each = 3)
  g <- trchlthgap(lifeexp = le)
  expect_identical(unique(g$period), 2021)
  b <- trchlthgap(lifeexp = le, infmort = .iim())
  expect_identical(unique(b$period), c("2021", "2014/2016"))
  le$estimate <- NA
  expect_error(trchlthgap(lifeexp = le), "no table given measures")
})
