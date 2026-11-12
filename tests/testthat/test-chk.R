# The argument checks now live in the compiled binding (src/rchk.cpp); each
# one is exercised through an exported function that uses it.

test_that("the number check names the argument and the rule broken", {
  expect_identical(incap(c(1L, 2L), 0.5, 1)[["lam"]], c(1, 2))
  expect_error(losbnd(numeric(0), 1), "`done` must be finite")
  expect_error(losbnd(c(1, Inf), 1), "`done` must be finite")
  expect_error(losbnd("1", 1), "`done` must be finite")
  expect_error(rtmhs(c(1, 2, 3), c(1, 2), 0), "length 3, not 2")
  expect_error(losbnd(-1, 1), "at least 0")
  expect_error(cllr(0, 1), "greater than 0")
  expect_identical(incap(1, 0.5, 0)[["s"]], 0)
})

test_that("the one-number check takes one number strictly inside", {
  expect_identical(hzdsel(0.5, 1, 0, c(1, 1), c(1, 1))[["share0"]], 0.5)
  expect_error(orrisk(2, 0, 0.5), "strictly between 0 and 1")
  expect_error(orrisk(2, c(0.2, 0.3), 0.5), "`rate`")
  expect_error(orrisk(2, NA_real_, 0.5), "`rate`")
})

test_that("the reference finds a group by name or by position", {
  expect_identical(expbnd(c(a = 1, b = 2), c(1, 1), 1, ref = "b")[["ratio"]],
                   c(0.5, 1))
  expect_identical(expbnd(c(a = 1, b = 2), c(1, 1), 1, ref = 2)[["ratio"]],
                   c(0.5, 1))
  expect_error(vodtst(c(1, 1), c(1, 1), ref = 1.5), "position from 1 to 2")
  expect_error(vodtst(c(1, 1), c(1, 1), ref = NA_character_), "`ref`")
  expect_identical(vodtst(c(4, 5), c(1, 1))[["group"]], c("1", "2"))
  expect_identical(vodtst(c(a = 4, b = 5), c(1, 1))[["group"]], c("a", "b"))
})
