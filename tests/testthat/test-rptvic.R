# Repeat victimisation. The kernel is checked by hand in src/t_rpv.cpp;
# here the binding against table(), dpois() and diff() in base R.

test_that("measures follow the counts by target", {
  set.seed(8)
  hh <- sample(1:400, 300, TRUE, prob = rgamma(400, 0.5))
  day <- sample(1:365, 300, TRUE)
  r <- rptvic(hh, t = day, pop = 1000, brk = c(0, 7, 30, 90, 366))
  expect_s3_class(r, "rptvic")
  k <- table(hh)
  v <- table(factor(k, 1:max(k)))
  expect_equal(r$tab$victims, c(1000 - length(k), as.numeric(v)))
  expect_equal(r$tab$incidents, r$tab$k * r$tab$victims)
  expect_equal(
    r$tab$expected,
    1000 * dpois(0:max(k), 300 / 1000),
    tolerance = 1e-13
  )
  s <- r$stat
  expect_equal(s[["victims"]], length(k))
  expect_equal(s[["conc"]], 300 / length(k))
  expect_equal(s[["rvict"]], sum(k >= 2))
  expect_equal(s[["erep"]], 1000 * (1 - ppois(1, 0.3)), tolerance = 1e-12)
  g <- unlist(lapply(split(day, hh), function(x) diff(sort(x))))
  expect_equal(r$gap[["median"]], median(g))
  expect_equal(r$gap[["mean"]], mean(g), tolerance = 1e-14)
  expect_equal(
    r$lag$n,
    as.numeric(table(cut(g, c(0, 7, 30, 90, 366), right = FALSE)))
  )
  expect_equal(r$lag$cum[4], 1)
  expect_output(print(r), "300 incidents")
})

test_that("without pop and t the rates per target and the gaps are absent", {
  r <- rptvic(c("a", "b", "a"))
  expect_identical(r$tab$k, 1:2)
  expect_true(is.na(r$stat[["prev"]]))
  expect_null(r$gap)
  expect_null(r$lag)
  expect_null(r$tab$expected)
})

test_that("rptvic() refuses bad input in words", {
  errs(
    "rptvic",
    list(
      list(
        "`pop` must be at least the number of victims",
        c("a", "b"),
        pop = 1
      ),
      list("`brk` needs the times", c("a", "b"), brk = c(0, 1)),
      list("`brk` must increase", c("a", "a"), t = c(1, 2), brk = c(2, 1)),
      list("`t`", c("a", "a"), t = c(1, NA))
    )
  )
})
