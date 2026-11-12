test_that(".num names the argument and the rule broken", {
  expect_identical(.num(c(1L, 2L), "x"), c(1, 2))
  expect_error(.num(numeric(0), "x"), "`x` must be finite")
  expect_error(.num(c(1, Inf), "x"), "`x` must be finite")
  expect_error(.num("1", "x"), "`x` must be finite")
  expect_error(.num(c(1, 2), "x", len = 3L), "length 3, not 2")
  expect_error(.num(-1, "x", lo = 0), "at least 0")
  expect_error(.num(0, "x", lo = 0, open = TRUE), "greater than 0")
  expect_identical(.num(0, "x", lo = 0), 0)
})

test_that(".one takes one number strictly inside the bounds", {
  expect_identical(.one(0.5, "r", 0, 1), 0.5)
  expect_error(.one(0, "r", 0, 1), "strictly between 0 and 1")
  expect_error(.one(c(0.2, 0.3), "r", 0, 1), "`r`")
  expect_error(.one(NA_real_, "r", 0, 1), "`r`")
})

test_that(".ref finds a group by name or by position", {
  expect_identical(.ref("b", c("a", "b"), 2L), 2L)
  expect_identical(.ref(2, c("a", "b"), 2L), 2L)
  expect_error(.ref(1.5, c("a", "b"), 2L), "position from 1 to 2")
  expect_error(.ref(NA_character_, c("a", "b"), 2L), "`ref`")
  expect_identical(.grp(c(4, 5)), c("1", "2"))
  expect_identical(.grp(c(a = 4, b = 5)), c("a", "b"))
})
