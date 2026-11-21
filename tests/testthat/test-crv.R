# cramerv() (R/crv.R): the association of a counted table, against
# ctabx() and against the formulas written out.

tab <- matrix(
  c(30, 6, 2, 9, 14, 3, 1, 4, 11),
  3,
  dimnames = list(c("none", "knife", "gun"), c("hands", "CEW", "firearm"))
)

test_that("cramerv() is ctabx()'s V, corrected or not", {
  s <- ctabx(tab)[["stat"]]
  r <- cramerv(tab)
  expect_s3_class(r, "cramerv")
  expect_named(
    r,
    c("v", "vraw", "lo", "hi", "chisq", "df", "p", "n", "minexp", "drop")
  )
  expect_identical(r[["v"]], s[["vbc"]])
  expect_identical(r[["vraw"]], s[["v"]])
  expect_identical(cramerv(tab, biascor = FALSE)[["v"]], s[["v"]])
  # the formula: sqrt(chisq / (n (k - 1)))
  ch <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  expect_equal(r[["chisq"]], unname(ch$statistic), tolerance = 1e-12)
  expect_equal(
    r[["vraw"]],
    sqrt(unname(ch$statistic) / (sum(tab) * 2)),
    tolerance = 1e-12
  )
  expect_equal(r[["minexp"]], min(ch$expected), tolerance = 1e-12)
  expect_true(r[["lo"]] <= r[["vraw"]] && r[["vraw"]] <= r[["hi"]])
  # a table() is taken as a matrix is
  expect_identical(unclass(cramerv(as.table(tab)))[1:8], unclass(r)[1:8])
  # independence: the corrected V is 0
  expect_identical(cramerv(matrix(c(10, 20, 10, 20), 2))[["v"]], 0)
  # a wider interval at a higher level
  w <- cramerv(tab, level = 0.99)
  expect_true(w[["lo"]] <= r[["lo"]] && w[["hi"]] >= r[["hi"]])
})

test_that("cramerv() drops empty margins and says when nothing is left", {
  z <- rbind(cbind(tab, none = 0), other = 0)
  r <- cramerv(z)
  expect_identical(r[["drop"]], 2)
  expect_equal(r[["vraw"]], cramerv(tab)[["vraw"]], tolerance = 1e-12)
  one <- cramerv(matrix(c(5, 7), 1))
  expect_true(is.na(one[["v"]]))
  expect_identical(one[["minexp"]], 5)
})

test_that("cramerv() prints its reading", {
  expect_output(
    print(cramerv(tab)),
    "bias-corrected; uncorrected",
    fixed = TRUE
  )
  expect_output(print(cramerv(tab)), "report the exact test", fixed = TRUE)
  expect_output(
    print(cramerv(tab, biascor = FALSE)),
    "(uncorrected; 95% interval",
    fixed = TRUE
  )
  big <- matrix(c(300, 60, 90, 140), 2)
  expect_false(any(grepl("exact test", capture.output(print(cramerv(big))))))
  expect_output(
    print(cramerv(rbind(tab, 0))),
    "1 empty rows or columns left out",
    fixed = TRUE
  )
  expect_invisible(print(cramerv(big)))
})

test_that("cramerv() refuses what is not a counted two-way table", {
  expect_error(
    cramerv(1:4),
    "cramerv: `tbl` must be a two-way table",
    fixed = TRUE
  )
  expect_error(cramerv(matrix(letters[1:4], 2)), "cramerv: `tbl`", fixed = TRUE)
  expect_error(cramerv(array(1, c(2, 2, 2))), "cramerv: `tbl`", fixed = TRUE)
  expect_error(
    cramerv(tab, biascor = NA),
    "cramerv: `biascor` must be TRUE or FALSE",
    fixed = TRUE
  )
  # a refusal of ctabx() comes under cramerv()'s name
  expect_error(cramerv(tab, level = 2), "^cramerv: ")
  expect_error(cramerv(matrix(c(1, -2, 3, 4), 2)), "^cramerv: ")
})
