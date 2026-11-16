test_that("tabrpt() gives the check, change, rates and trend of a table", {
  d <- data.frame(year = rep(2019:2023, 2),
                  svc = rep(c("A", "B"), each = 5),
                  n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33),
                  pop = rep(c(50000, 12000), each = 5))
  r <- tabrpt(d, "n", "year", by = "svc", pop = "pop", rnd = 5)
  expect_s3_class(r, "tabrpt")
  expect_identical(names(r), c("check", "change", "rates", "trend"))
  # each step is the function it calls, on the columns written out
  expect_identical(r$check, catchk(d["svc"]))
  expect_identical(r$change, yoypx(d$n, d$year, by = d$svc, rnd = 5))
  expect_identical(r$rates, grprr(d$n, d$pop, grp = d$svc, per = 1000))
  expect_identical(r$trend$A, trndrr(d$n[1:5], 2019:2023,
                                     expo = rep(50000, 5)))
  expect_identical(r$trend$B, trndrr(d$n[6:10], 2019:2023,
                                     expo = rep(12000, 5)))
  expect_output(print(r), "^tabrpt: 10 group-periods")
})

test_that("tabrpt() sums a period's rows and passes over short groups", {
  d <- data.frame(t = c(1, 1, 2, 3, 4, 1, 2),
                  g1 = c("a", "a", "a", "a", "a", "b", "b"),
                  g2 = "x", n = c(3, 4, 6, 9, 12, 5, 6))
  r <- tabrpt(d, "n", "t", by = c("g1", "g2"))
  expect_identical(names(r$trend), c("a / x", "b / x"))
  # written out: the two rows of period 1 of group a summed, 3 + 4 = 7
  expect_identical(r$trend[["a / x"]], trndrr(c(7, 6, 9, 12), c(1, 2, 3, 4)))
  expect_null(r$trend[["b / x"]])
  expect_null(r$rates)
  one <- tabrpt(d[1:5, ], "n", "t")
  expect_null(one$check)
  expect_identical(names(one$trend), "all")
})

test_that("tabrpt() refuses what it cannot read, and says so as itself", {
  d <- data.frame(t = 1:3, n = c(1, 2, 3))
  expect_error(tabrpt(list(t = 1), "n", "t"),
               "tabrpt: `data` must be a data frame with rows.", fixed = TRUE)
  expect_error(tabrpt(d, "m", "t"),
               "tabrpt: `value` must name one column of `data`.",
               fixed = TRUE)
  expect_error(tabrpt(d, "n", c("t", "n")),
               "tabrpt: `prd` must name one column of `data`.", fixed = TRUE)
  expect_error(tabrpt(d, "n", "t", by = "z"),
               "tabrpt: `by` must name columns of `data`.", fixed = TRUE)
  expect_error(tabrpt(d, "n", "t", pop = 3),
               "tabrpt: `pop` must name one column of `data`.", fixed = TRUE)
  expect_error(tabrpt(transform(d, n = c(1, -2, 3)), "n", "t"), "^tabrpt: ")
})
