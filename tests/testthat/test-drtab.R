# drfttab() and drftcal() (R/drtab.R): the column tests are drftks(),
# drftpsi() and drfthomog(), checked here to be the ones reported.

set.seed(11)
ref <- data.frame(
  age = stats::rnorm(300, 35, 9),
  charge = sample(c("assault", "theft", "drugs"), 300, TRUE),
  stringsAsFactors = FALSE
)
cur <- data.frame(
  age = stats::rnorm(300, 40, 9),
  charge = sample(c("assault", "theft", "drugs"), 300, TRUE),
  extra = 1,
  stringsAsFactors = FALSE
)

test_that("drfttab() reports the column tests of the family", {
  r <- drfttab(ref, cur)
  expect_s3_class(r, c("drfttab", "data.frame"))
  expect_named(
    r,
    c(
      "column",
      "type",
      "test",
      "statistic",
      "p",
      "psi",
      "nref",
      "ncur",
      "flag",
      "note"
    )
  )
  expect_identical(r[["column"]], c("age", "charge"))
  expect_identical(r[["type"]], c("numeric", "label"))
  k <- drftks(ref$age, cur$age)
  expect_identical(r[["statistic"]][1L], k[["statistic"]])
  expect_identical(r[["p"]][1L], k[["p_value"]])
  expect_identical(r[["psi"]][1L], drftpsi(ref$age, cur$age)[["psi"]])
  h <- drfthomog(ref$charge, cur$charge)
  expect_identical(r[["statistic"]][2L], h[["statistic"]])
  expect_identical(r[["p"]][2L], h[["p_value"]])
  expect_true(is.na(r[["psi"]][2L]))
  expect_identical(r[["nref"]], c(300L, 300L))
  expect_true(r[["flag"]][1L])
  expect_identical(attr(r, "only"), list(ref = character(), cur = "extra"))
  expect_identical(attr(r, "alpha"), 0.01)
  # the bins reach drftpsi()
  expect_identical(
    drfttab(ref, cur, bins = 4L)[["psi"]][1L],
    drftpsi(ref$age, cur$age, bins = 4L)[["psi"]]
  )
})

test_that("drfttab() flags by the test or the index, not a small sample's", {
  a <- data.frame(x = seq(0, 1, length.out = 120))
  b <- data.frame(x = seq(0, 1, length.out = 120) + 0.22)
  r <- drfttab(a, b, alpha = 1e-12)
  expect_gt(r[["p"]], 1e-12)
  expect_gt(r[["psi"]], 0.25)
  expect_true(r[["flag"]])
  expect_false(drfttab(a, b, alpha = 1e-12, psiminn = 500L)[["flag"]])
  expect_false(drfttab(a, b, alpha = 1e-12, psithr = 50)[["flag"]])
  # unchanged data is not flagged
  expect_false(drfttab(a, a)[["flag"]])
})

test_that("drfttab() says why a column was not tested", {
  a <- data.frame(
    n = c(NA_real_, NA),
    l = c(NA_character_, NA),
    d = as.Date(c("2026-01-01", "2026-03-01")),
    one = c("a", "a"),
    mix = c(1, 2),
    stringsAsFactors = FALSE
  )
  a$odd <- list(1, 2)
  b <- data.frame(
    n = c(1, 2),
    l = c("a", "b"),
    d = as.Date(c("2026-02-01", "2026-04-01")),
    one = c("a", "a"),
    mix = c("x", "y"),
    stringsAsFactors = FALSE
  )
  b$odd <- list(3, 4)
  r <- drfttab(a, b)
  n <- stats::setNames(r[["note"]], r[["column"]])
  f <- stats::setNames(r[["flag"]], r[["column"]])
  expect_identical(n[["n"]], "no values on one side")
  expect_identical(n[["l"]], "no values on one side")
  expect_true(is.na(f[["n"]]) && is.na(f[["l"]]))
  # a date column is tested as numbers
  expect_identical(r[["type"]][r[["column"]] == "d"], "numeric")
  expect_false(is.na(r[["p"]][r[["column"]] == "d"]))
  # one category: the homogeneity test does not apply, and says so
  expect_match(n[["one"]], "inapplicable", fixed = TRUE)
  expect_true(is.na(f[["one"]]))
  expect_identical(
    n[["mix"]],
    "not numbers on both sides nor labels on both sides"
  )
  expect_identical(r[["type"]][r[["column"]] == "odd"], "other")
  expect_output(print(r), "5 not testable", fixed = TRUE)
})

test_that("drfttab() prints, with the columns of one table alone", {
  expect_output(
    print(drfttab(ref, cur)),
    "drfttab: 2 columns, 1 flagged, 0 not testable",
    fixed = TRUE
  )
  expect_output(print(drfttab(ref, cur)), "only in `cur`: extra", fixed = TRUE)
  expect_output(print(drfttab(cur, ref)), "only in `ref`: extra", fixed = TRUE)
  expect_invisible(print(drfttab(ref, ref)))
})

test_that("drfttab() refuses bad input", {
  expect_error(
    drfttab(1, cur),
    "drfttab: `ref` and `cur` must be data frames",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, "x"),
    "drfttab: `ref` and `cur` must be data frames",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, data.frame(other = 1)),
    "drfttab: `ref` and `cur` share no column",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, cur, alpha = 1),
    "drfttab: `alpha` must be one number in (0, 1)",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, cur, alpha = c(0.1, 0.2)),
    "drfttab: `alpha`",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, cur, psithr = 0),
    "drfttab: `psithr` must be one positive number",
    fixed = TRUE
  )
  expect_error(
    drfttab(ref, cur, psiminn = 1.5),
    "drfttab: `psiminn` must be a whole number",
    fixed = TRUE
  )
  # a refusal of a column test comes under drfttab()'s name
  expect_error(drfttab(ref, cur, bins = 1), "^drfttab: ")
})

test_that("drftcal() measures the false alarms of unchanged data", {
  before <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  r <- drftcal(ref, n = 12, alpha = 0.2)
  after <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  expect_identical(before, after)
  expect_s3_class(r, "drftcal")
  expect_named(r, c("columns", "anyflag", "alpha", "n", "alphafam"))
  expect_identical(r[["columns"]][["column"]], c("age", "charge"))
  expect_true(all(r[["columns"]][["far"]] >= 0 & r[["columns"]][["far"]] <= 1))
  expect_identical(r[["n"]], 12L)
  expect_identical(r[["alphafam"]], 0.1)
  expect_true(r[["anyflag"]] >= max(r[["columns"]][["far"]]))
  # the splits are those of the kernel's generator: split 1 by hand
  u <- cunif(nrow(ref), floor(cunif(12, 1)[1L] * 2^52))
  k <- sort(order(u)[seq_len(nrow(ref) %/% 2L)])
  one <- drfttab(ref[k, ], ref[-k, ], alpha = 0.2)
  expect_identical(
    drftcal(ref, n = 1, alpha = 0.2)[["columns"]][["far"]],
    as.numeric(one[["flag"]])
  )
  # the same seed, the same splits
  expect_identical(drftcal(ref, n = 12, alpha = 0.2), r)
  expect_false(identical(drftcal(ref, n = 12, alpha = 0.2, seed = 2L), r))
  # the further arguments reach drfttab()
  s <- drftcal(ref, n = 4, alpha = 1e-12, psithr = 1e-9, psiminn = 1L)
  expect_identical(s[["columns"]][["far"]][1L], 1)
  expect_output(print(s), "any column flagged on 100% of splits", fixed = TRUE)
  expect_output(print(s), "age", fixed = TRUE)
  q <- drftcal(ref, n = 3, alpha = 1e-12)
  expect_identical(q[["anyflag"]], 0)
  expect_output(print(q), "no column was flagged on any split", fixed = TRUE)
  expect_invisible(print(q))
})

test_that("drftcal() refuses bad input", {
  expect_error(
    drftcal(1),
    "drftcal: `data` must be a data frame with at least four rows",
    fixed = TRUE
  )
  expect_error(
    drftcal(ref[1:3, ]),
    "drftcal: `data` must be a data frame",
    fixed = TRUE
  )
  expect_error(
    drftcal(ref, n = 0),
    "drftcal: `n` must be a whole number from 1 to 1000000",
    fixed = TRUE
  )
  expect_error(drftcal(ref, n = 2.5), "drftcal: `n`", fixed = TRUE)
  expect_error(drftcal(ref, seed = "a"), "drftcal: `seed`", fixed = TRUE)
  expect_error(
    drftcal(ref, seed = -1),
    "drftcal: `seed` must be one whole number from 0, below 2^53",
    fixed = TRUE
  )
  expect_error(drftcal(ref, n = 1e7), "drftcal: `n`", fixed = TRUE)
  expect_error(drftcal(ref, n = 2, alpha = 2), "^drftcal: ")
})
