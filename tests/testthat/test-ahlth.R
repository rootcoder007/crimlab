# Air pollution and health (R/ahlth.R): the band a baseline rate per
# person-year must fall in, and the burden of each area of a table.

test_that("a baseline rate outside the band is refused unless declared", {
  expect_error(mortdisp(1, 1000, 800, 0.01), "looks like a rate per 100,000")
  expect_error(mortdisp(1, 1000, 1e-12, 0.01),
               "below one case per 100 billion")
  expect_error(mortdisp(1, 1000, 1e-6, 0.01), "rare = TRUE")
  expect_error(mortdisp(1, 1000, 0.2, 0.01), "hirate = TRUE")
  expect_error(mortdisp(1, 1000, 0.008, 0.01, rare = NA),
               "mortdisp: `rare` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(mortdisp(1, 1000, 0.008, 0.01, hirate = 1),
               "`hirate` must be TRUE or FALSE")
  expect_identical(mortdisp(1, 1000, 1e-6, 0.01, rare = TRUE),
                   mortdisp(1, 1000, 1e-6, 0.01, rare = TRUE, hirate = TRUE))
  expect_identical(mortdisp(1, 1000, 0, 0.01), mortdisp(2, 1000, 0, 0.01))
})

test_that("the burden of each area, worst first, with its columns' names", {
  t <- data.frame(fsa = c("A1", "B2"), exposure = c(18, 28),
                  population = c(60000, 40000), brate = 0.008)
  r <- polbrdarea(t)
  expect_s3_class(r, "polbrdarea")
  expect_identical(r$fsa, c("B2", "A1"))
  expect_named(r, c("fsa", "exposure", "population", "brate", "rr", "paf",
                    "attrcases", "basecases"))
  b <- polbrdn(28, 1, 0.008, 40000)
  expect_identical(r$attrcases[1L], b$attrcases)
  t$population[1L] <- -1
  expect_error(polbrdarea(t), "`popcol`")
  t$population <- c("a", "b")
  expect_error(polbrdarea(t),
               "column `population` must be numeric (a factor or character",
               fixed = TRUE)
  expect_error(polbrdarea(t[0L, ]), "areatab has no rows")
  expect_error(polbrdarea(t[c("fsa", "exposure")]),
               "areatab missing columns: population, brate", fixed = TRUE)
})
