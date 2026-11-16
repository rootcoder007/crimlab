test_that("pscred() credits each day at the ratio and leaves the rest", {
  d <- c(10, 30, 45)
  r <- pscred(d)
  expect_identical(names(r), c("days", "ratio", "credit"))
  expect_identical(r$credit, d * 1.5)
  s <- pscred(d, ratio = c(1, 1.5, 1.5), sentence = c(60, 40, 60))
  # written out: 10, 45, 67.5 days of credit; 60 - 10, 40 - 45 (none
  # left), 60 - 67.5 (none left)
  expect_identical(s$credit, c(10, 45, 67.5))
  expect_identical(s$remain, c(50, 0, 0))
  expect_identical(s$sentence, c(60, 40, 60))
  expect_identical(pscred(0, ratio = 0, sentence = 5)$remain, 5)
})

test_that("pscred() refuses what is not a count of days or a ratio", {
  expect_error(pscred(-1), "pscred: `days` must be at least 0.",
               fixed = TRUE)
  expect_error(pscred("a"), "pscred: `days` must be numeric, got character.",
               fixed = TRUE)
  expect_error(pscred(NA_real_),
               "pscred: `days` must be finite, with no missing values.",
               fixed = TRUE)
  expect_error(pscred(numeric(0)),
               "pscred: `days` must have at least one value.", fixed = TRUE)
  expect_error(pscred(10, ratio = 2),
               "pscred: `ratio` must be between 0 and 1.5.", fixed = TRUE)
  expect_error(pscred(c(1, 2, 3), ratio = c(1, 1.5)),
               paste("pscred: `ratio` must have one value or one per value",
                     "of `days`."), fixed = TRUE)
  expect_error(pscred(1, sentence = -2),
               "pscred: `sentence` must be at least 0.", fixed = TRUE)
  expect_error(pscred(c(1, 2, 3), sentence = c(1, 2)),
               paste("pscred: `sentence` must have one value or one per",
                     "value of `days`."), fixed = TRUE)
})
