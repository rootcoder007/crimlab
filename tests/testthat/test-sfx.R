test_that("synfix() draws a preset as sample() draws after set.seed()", {
  x <- synfix("otis-b02", n = 20, seed = 3)
  # written out: the columns of the published OTIS b02 table, each drawn
  # with sample() in turn
  set.seed(3)
  want <- data.frame(
    EndFiscalYear = sample(2018:2024, 20, replace = TRUE),
    UniqueIndividual_ID = sprintf("syn-%05d",
                                  sample.int(5, 20, replace = TRUE)),
    Gender = sample(c("Female", "Male"), 20, replace = TRUE),
    Region_MostRecentPlacement = sample(
      c("Central", "Eastern", "Northern", "Toronto", "Western"), 20,
      replace = TRUE
    ),
    Age_Category = sample(c("18 to 24", "25 to 49", "50+"), 20,
                          replace = TRUE),
    TotalAggregatedDays_Segregation = sample(0:80, 20, replace = TRUE),
    stringsAsFactors = FALSE
  )
  expect_identical(x, want)
})

test_that("synfix() gives a family as a list, one seed after another", {
  x <- synfix("cuof", n = 5, seed = 2)
  expect_length(x, 12L)
  expect_identical(names(x)[1:3], c("incidents", "incident_type",
                                    "indigenous"))
  set.seed(3)
  want <- data.frame(
    year = (2023:2024)[sample.int(2, 5, replace = TRUE)],
    incident_ID = sprintf("syn-%05d", sample.int(2, 5, replace = TRUE)),
    status = "Final"[sample.int(1, 5, replace = TRUE)],
    stringsAsFactors = FALSE
  )
  expect_identical(x[["incidents"]], want)
  expect_identical(x[["race"]], synfix("cuof-race", n = 5, seed = 2 + 8))
})

test_that("synfix() draws a table given as lists of values", {
  x <- synfix(list(a = c("x", "y", "z"), b = 5:9, c = 1.5), n = 10,
              seed = 7)
  set.seed(7)
  want <- data.frame(a = c("x", "y", "z")[sample.int(3, 10, TRUE)],
                     b = (5:9)[sample.int(5, 10, TRUE)],
                     c = 1.5[sample.int(1, 10, TRUE)])
  expect_identical(x, want)
  p <- synfix()
  expect_named(p, c("kind", "column", "type", "values"))
  expect_length(unique(p$kind), 41L)
  expect_identical(p$values[p$kind == "otis-b02" & p$column == "Gender"],
                   "Female|Male")
})

test_that("synfix() refuses what it cannot draw", {
  expect_error(synfix(3),
               paste("synfix: `kind` must be the name of a preset, a named",
                     "list of the values of each column, or NULL."),
               fixed = TRUE)
  expect_error(synfix("otis-z99"),
               paste("synfix: `kind` names no preset: \"otis-z99\";",
                     "synfix() lists them."), fixed = TRUE)
  expect_error(synfix(list(1:3)),
               "synfix: a list `kind` must name each column once.",
               fixed = TRUE)
  expect_error(synfix(list(a = c(1, NA))),
               paste("synfix: column `a` of `kind` must be a vector of the",
                     "values it takes, with no missing values."),
               fixed = TRUE)
  expect_error(synfix("otis-a01", n = 0),
               "synfix: `n` must be between 1 and 1e+07.", fixed = TRUE)
  expect_error(synfix("otis-a01", seed = -1),
               "synfix: `seed` must be between 0 and 2147483647.",
               fixed = TRUE)
  expect_error(synfix("otis-a01", n = 2.5),
               "synfix: `n` must be a whole number.", fixed = TRUE)
})
