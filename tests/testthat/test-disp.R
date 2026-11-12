test_that("vodtst matches the odds ratio by hand", {
  r <- vodtst(c(w = 60, b = 40), c(w = 50, b = 50))
  expect_equal(r[["ratio"]], c(1, (40 / 60) / (50 / 50)))
  expect_equal(r[["share_day"]], c(0.6, 0.4))
  expect_equal(r[["sign"]], c(1L, -1L))
  expect_false("adj" %in% names(r))
  # an exposure shift equal to the ratio explains all of it
  a <- vodtst(c(w = 60, b = 40), c(w = 50, b = 50), shift = c(1, 2 / 3))
  expect_equal(a[["adj"]], c(1, 1))
  # the reference by name moves the baseline
  expect_equal(vodtst(c(w = 60, b = 40), c(w = 50, b = 50), ref = "b")[["ratio"]],
               c(1.5, 1))
  expect_equal(vodtst(c(5, 5), c(7, 7))[["sign"]], c(0L, 0L))
})

test_that("orrisk returns the pair with that odds ratio and that rate", {
  for (or in c(0.05, 0.5, 1, 3, 400)) {
    for (q in c(0.01, 0.4, 0.99)) {
      for (s in c(0.02, 0.5, 0.9)) {
        r <- orrisk(or, q, s)
        a <- r[["exposed"]]
        b <- r[["unexposed"]]
        expect_true(a > 0 && a < 1 && b > 0 && b < 1)
        expect_equal(a * (1 - b) / (b * (1 - a)), or, tolerance = 1e-9)
        expect_equal(s * a + (1 - s) * b, q, tolerance = 1e-12)
        # Csls.RelRisk.or_overstates: the relative risk sits between 1 and OR
        expect_true(abs(log(r[["rr"]])) <= abs(log(or)) + 1e-12)
      }
    }
  }
})

test_that("outtst shows infra-marginality", {
  expect_equal(outtst(c(0.1, 0.5, 0.9), 0.3), c(stop = 2 / 3, hit = 0.7))
  expect_equal(outtst(c(0.1, 0.5, 0.9), 0.3, w = c(1, 3, 1)),
               c(stop = 0.8, hit = (1.5 + 0.9) / 4))
  expect_true(is.nan(outtst(c(0.1, 0.2), 0.5)[["hit"]]))
  # same threshold, different hit rates
  expect_gt(outtst(c(0.2, 0.9), 0.5)[["hit"]], outtst(c(0.2, 0.6), 0.5)[["hit"]])
})

test_that("dspbm is a product", {
  r <- dspbm(c(a = 9000, b = 1000), c(a = 900, b = 300), c(a = 9, b = 6))
  expect_equal(r[["cd"]], c(1, 3))
  expect_equal(r[["fd"]], c(1, 2))
  expect_equal(r[["rd"]], r[["cd"]] * r[["fd"]])
})

test_that("dyadrt is flat under random mixing", {
  pop <- c(a = 800, b = 150, c = 50)
  r <- dyadrt(outer(pop, pop) * 0.002, pop)
  expect_equal(r[["ratio"]], rep(1, 9))
  expect_equal(r[["rate_off"]], r[["null_off"]])
  # per offender rises with the victim group, per pair does not
  expect_gt(r[["rate_off"]][1], r[["rate_off"]][7])
  expect_equal(r[["rate_pair"]][1], r[["rate_pair"]][7])
  expect_equal(r[["off"]][2], "b")
  expect_equal(r[["vic"]][2], "a")
})

test_that("bad input stops in R with the argument named", {
  expect_error(vodtst(c(1, NA), c(1, 1)), "`day`")
  expect_error(vodtst(c(1, 1), c(1, 0)), "`night`")
  expect_error(vodtst(c(1, 1), c(1, 1, 1)), "length 2")
  expect_error(vodtst(c(0, 1), c(1, 1)), "no daylight")
  expect_error(vodtst(c(1, 1), c(1, 1), ref = 3), "`ref`")
  expect_error(vodtst(c(a = 1, b = 1), c(1, 1), ref = "z"), "names no group")
  expect_error(vodtst(c(1, 1), c(1, 1), shift = c(1, -1)), "`shift`")
  expect_error(orrisk(0, 0.5, 0.5), "`or`")
  expect_error(orrisk(2, 1, 0.5), "`rate`")
  expect_error(orrisk(2, 0.5, 0), "`share`")
  expect_error(outtst(c(0.5, 1.5), 0.2), "at most 1")
  expect_error(outtst(c(0.5, 0.6), 0.2, w = c(0, 0)), "positive total")
  expect_error(outtst("a", 0.2), "`q`")
  expect_error(dspbm(c(1, 1), c(1, 0), c(1, 1)), "`contact`")
  expect_error(dyadrt(matrix(1, 2, 3), c(1, 1)), "2 by 2")
  expect_error(dyadrt(matrix(0, 2, 2), c(1, 1)), "no offences")
  expect_error(dyadrt(matrix(-1, 2, 2), c(1, 1)), "`count`")
})

test_that("proof names real theorems", {
  p <- proof()
  expect_true(all(c("vodtst", "orrisk", "outtst", "dspbm", "dyadrt", "gtmp",
                    "kinidx") %in% p[["fn"]]))
  expect_equal(nrow(proof("dspbm")), 2L)
  expect_error(proof("nope"), "no proof recorded")
  expect_error(proof(1), "function names")
  # every function proof() names is exported
  expect_true(all(p[["fn"]] %in% getNamespaceExports("csls")))
  for (k in seq_len(nrow(p))) {
    f <- system.file(p[["file"]][k], package = "csls")
    expect_true(nzchar(f), info = p[["file"]][k])
    src <- readLines(f, warn = FALSE)
    nm <- sub(".*\\.", "", p[["theorem"]][k])
    expect_true(any(grepl(paste0("^theorem ", nm, "\\b"), src)), info = nm)
    # the namespace in the file is the one proof() gives
    ns <- sub("\\.[^.]*$", "", p[["theorem"]][k])
    expect_true(any(src == paste("namespace", ns)), info = ns)
  }
  expect_false(any(grepl("sorry", unlist(lapply(
    list.files(system.file("lean", package = "csls"), full.names = TRUE),
    readLines, warn = FALSE
  )))))
})

test_that("results agree with rmorie", {
  skip_if_not_installed("rmorie")
  day <- c(w = 600, b = 400)
  night <- c(w = 700, b = 300)
  m <- rmorie::morie_veil_of_darkness(day, night)
  v <- vodtst(day, night)
  expect_equal(v[["ratio"]], m[["bias_ratio"]], tolerance = 1e-14)
  expect_equal(v[["share_day"]], m[["share_day"]], tolerance = 1e-14)

  m <- rmorie::morie_relative_risk_from_or(3, 0.4, 0.5)[["risks"]]
  r <- orrisk(3, 0.4, 0.5)
  # rmorie finds the root numerically; the closed form here is the reference
  expect_equal(unname(r[c("exposed", "unexposed", "rr")]), unname(m),
               tolerance = 1e-6)

  g <- list(A = list(q = c(0.3, 0.9), w = c(1, 1), t = 0.5),
            B = list(q = 0.9, w = 1, t = 0.2))
  m <- rmorie::morie_outcome_test(g)
  for (k in seq_along(g)) {
    o <- outtst(g[[k]][["q"]], g[[k]][["t"]], g[[k]][["w"]])
    expect_equal(o[["stop"]], m[["stop_rate"]][k], tolerance = 1e-14)
    expect_equal(o[["hit"]], m[["hit_rate"]][k], tolerance = 1e-14)
  }

  pop <- c(a = 9000, b = 1000)
  con <- c(a = 900, b = 300)
  frc <- c(a = 9, b = 6)
  m <- rmorie::morie_disparity_benchmark(pop, con, frc, reference = "a")
  d <- dspbm(pop, con, frc)
  expect_equal(d[["cd"]], m[["contact_disparity"]], tolerance = 1e-14)
  expect_equal(d[["fd"]], m[["force_given_contact_disparity"]],
               tolerance = 1e-14)
  expect_equal(d[["rd"]], m[["resident_disparity"]], tolerance = 1e-14)

  pop <- c(A = 80000, B = 20000)
  m <- rmorie::morie_interracial_rates(
    c(A_on_A = 900, B_on_A = 200, A_on_B = 120, B_on_B = 300), pop
  )
  d <- dyadrt(matrix(c(900, 200, 120, 300), 2), pop)
  expect_equal(d[["off"]], m[["offender"]])
  expect_equal(d[["vic"]], m[["victim"]])
  expect_equal(d[["rate_off"]], m[["rate_per_offender_group"]],
               tolerance = 1e-13)
  expect_equal(d[["rate_pair"]], m[["rate_per_pair_exposure"]],
               tolerance = 1e-13)
  expect_equal(d[["null_off"]], m[["null_rate_per_offender_group"]],
               tolerance = 1e-13)
  expect_equal(d[["ratio"]], m[["ratio_to_null"]], tolerance = 1e-13)
})
