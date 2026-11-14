# Values of the disparity kernels are checked in C++ (src/t_disp.cpp); here
# what only R can show: labels, messages, the proofs and rmorie parity
# (rmorie's values stored in ref/disp.txt by dev/ref/disp.R, version named
# there).

test_that("vodtst and dyadrt label what they return", {
  r <- vodtst(c(w = 60, b = 40), c(w = 50, b = 50))
  expect_false("adj" %in% names(r))
  pop <- c(a = 800, b = 150, c = 50)
  r <- dyadrt(outer(pop, pop) * 0.002, pop)
  expect_identical(r[["off"]][2], "b")
  expect_identical(r[["vic"]][2], "a")
})

test_that("orrisk: the relative risk sits between 1 and the odds ratio", {
  # Crimlab.RelRisk.or_overstates; rr is formed in R
  for (or in c(0.05, 0.5, 1, 3, 400)) {
    for (q in c(0.01, 0.4, 0.99)) {
      for (s in c(0.02, 0.5, 0.9)) {
        r <- orrisk(or, q, s)
        expect_lte(abs(log(r[["rr"]])), abs(log(or)) + 1e-12)
      }
    }
  }
})

test_that("bad input stops in R with the argument named", {
  errs("vodtst", list(
    list("`day`", c(1, NA), c(1, 1)),
    list("`night`", c(1, 1), c(1, 0)),
    list("length 2", c(1, 1), c(1, 1, 1)),
    list("no daylight", c(0, 1), c(1, 1)),
    list("`ref`", c(1, 1), c(1, 1), ref = 3),
    list("names no group", c(a = 1, b = 1), c(1, 1), ref = "z"),
    list("`shift`", c(1, 1), c(1, 1), shift = c(1, -1))
  ))
  errs("orrisk", list(
    list("`or`", 0, 0.5, 0.5),
    list("`rate`", 2, 1, 0.5),
    list("`share`", 2, 0.5, 0)
  ))
  errs("outtst", list(
    list("at most 1", c(0.5, 1.5), 0.2),
    list("positive total", c(0.5, 0.6), 0.2, w = c(0, 0)),
    list("`q`", "a", 0.2)
  ))
  errs("dspbm", list(list("`contact`", c(1, 1), c(1, 0), c(1, 1))))
  errs("dyadrt", list(
    list("2 by 2", matrix(1, 2, 3), c(1, 1)),
    list("no offences", matrix(0, 2, 2), c(1, 1)),
    list("`count`", matrix(-1, 2, 2), c(1, 1))
  ))
})

test_that("proof names real theorems", {
  p <- proof()
  expect_true(all(c("vodtst", "orrisk", "outtst", "dspbm", "dyadrt", "gtmp",
                    "kinidx") %in% p[["fn"]]))
  expect_identical(nrow(proof("dspbm")), 2L)
  expect_error(proof("nope"), "no proof recorded")
  expect_error(proof(1), "function names")
  # every function proof() names is exported
  expect_true(all(p[["fn"]] %in% getNamespaceExports("crimlab")))
  for (k in seq_len(nrow(p))) {
    f <- system.file(p[["file"]][k], package = "crimlab")
    expect_true(nzchar(f), info = p[["file"]][k])
    src <- readLines(f, warn = FALSE)
    nm <- sub(".*\\.", "", p[["theorem"]][k])
    expect_true(any(grepl(paste0("^theorem ", nm, "\\b"), src)), info = nm)
    # the namespace in the file is the one proof() gives
    ns <- sub("\\.[^.]*$", "", p[["theorem"]][k])
    expect_true(any(src == paste("namespace", ns)), info = ns)
  }
  expect_false(any(grepl("sorry", unlist(lapply(
    list.files(system.file("lean", package = "crimlab"), "\\.lean$",
               full.names = TRUE),
    readLines, warn = FALSE
  )), fixed = TRUE)))
})

test_that("results agree with rmorie", {
  ref <- refget("disp")
  day <- c(w = 600, b = 400)
  night <- c(w = 700, b = 300)
  # from rmorie::morie_veil_of_darkness(day, night)
  v <- vodtst(day, night)
  expect_equal(v[["ratio"]], ref$vod$ratio, tolerance = 1e-14)
  expect_equal(v[["shday"]], ref$vod$shday, tolerance = 1e-14)

  # from rmorie::morie_relative_risk_from_or(3, 0.4, 0.5)[["risks"]]
  r <- orrisk(3, 0.4, 0.5)
  # rmorie finds the root numerically; the closed form here is the reference
  expect_equal(unname(r[c("exposed", "unexposed", "rr")]), ref$orr,
               tolerance = 1e-6)

  g <- list(A = list(q = c(0.3, 0.9), w = c(1, 1), t = 0.5),
            B = list(q = 0.9, w = 1, t = 0.2))
  # from rmorie::morie_outcome_test(g)
  m <- ref$out
  for (k in seq_along(g)) {
    o <- outtst(g[[k]][["q"]], g[[k]][["t"]], g[[k]][["w"]])
    expect_equal(o[["stop"]], m$stop[k], tolerance = 1e-14)
    expect_equal(o[["hit"]], m$hit[k], tolerance = 1e-14)
  }

  pop <- c(a = 9000, b = 1000)
  con <- c(a = 900, b = 300)
  frc <- c(a = 9, b = 6)
  # from rmorie::morie_disparity_benchmark(pop, con, frc, reference = "a")
  d <- dspbm(pop, con, frc)
  expect_equal(d[["cd"]], ref$bm$cd, tolerance = 1e-14)
  expect_equal(d[["fd"]], ref$bm$fd, tolerance = 1e-14)
  expect_equal(d[["rd"]], ref$bm$rd, tolerance = 1e-14)

  pop <- c(A = 80000, B = 20000)
  # rmorie::morie_interracial_rates(c(A_on_A = 900, B_on_A = 200,
  #   A_on_B = 120, B_on_B = 300), pop)
  m <- ref$dy
  d <- dyadrt(matrix(c(900, 200, 120, 300), 2), pop)
  expect_identical(d[["off"]], m$off)
  expect_identical(d[["vic"]], m$vic)
  expect_equal(d[["rtoff"]], m$rtoff, tolerance = 1e-13)
  expect_equal(d[["rtpair"]], m$rtpair, tolerance = 1e-13)
  expect_equal(d[["nulloff"]], m$nulloff, tolerance = 1e-13)
  expect_equal(d[["ratio"]], m$ratio, tolerance = 1e-13)
})
