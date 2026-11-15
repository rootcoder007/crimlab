# print, format and summary of the scaling (svfit) and count-dynamics
# (tsfit) results, and the one print method that every result whose
# format method gives its lines shares (R/svfmt.R).

test_that("a scaling result prints its lines and returns itself unseen", {
  f <- cmdsf(as.matrix(dist(1:5)), 1)
  expect_identical(format(f)[1L], "cmdsf:")
  expect_identical(utils::capture.output(print(f)), format(f))
  p <- expect_output(expect_invisible(print(f)), "cmdsf:")
  expect_identical(p, f)
  s <- summary(f)
  expect_identical(s$parts$part, names(f))
  expect_output(expect_invisible(print(s)), "parts")
})

test_that("a count-dynamics result is printed the same way", {
  f <- lvfit(c(410, 455, 430, 470, 520, 498, 515))
  expect_identical(utils::capture.output(print(f)), format(f))
  expect_output(expect_invisible(print(f)), "lvfit:")
  expect_s3_class(summary(f), "summary.svfit")
})

test_that("every class the shared print method serves is registered", {
  for (k in c("cnsprof", "profshr", "cnscwb", "trchlthgap", "mhcont",
              "cwbscr", "connassoc", "healchg", "commidx", "frmaud",
              "sctidy", "limfrom", "indgap", "exptrnd", "trcprog", "ordcmp",
              "stdprev", "welprof", "holbal", "dualcont", "svfit", "tsfit")) {
    expect_identical(utils::getS3method("print", k), .prfmt)
  }
})
