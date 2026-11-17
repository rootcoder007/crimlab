# The person-year table (otpyr()) on small made-up tables with the columns
# of the restrictive confinement release: every expected value is worked
# out by hand in the comments. Nothing here reaches the network (otpyr() by
# key reads a copy of the service in a temporary folder, helper-od.R).

a01 <- function() {
  data.frame(
    EndFiscalYear = c(2025, 2025, 2025, 2025, 2026, 2026, 2026),
    UniqueIndividual_ID = c("a", "a", "a", "b", "a", "c", "c"),
    Region_AtTimeOfPlacement = c(
      "West",
      "East",
      "West",
      "East",
      "North",
      "East",
      "East"
    ),
    Region_MostRecentPlacement = c(
      "West",
      "West",
      "West",
      "East",
      "East",
      "East",
      "East"
    ),
    Gender = c("Male", "Male", "Male", "Female", "Male", "Male", "Male"),
    Age_Category = "25 to 49",
    MentalHealth_Alert = c("Yes", "No", "No", "No", "No", "Yes", "Yes"),
    SuicideRisk_Alert = c("No", "No", "No", "No", "No", "Yes", "Yes"),
    SuicideWatch_Alert = c("No", "No", "No", "No", "No", "Yes", "Yes"),
    Number_Of_Placements = c(3, 1, 2, 5, 4, 1, 1)
  )
}

test_that("otpyr() builds the person-year table from the release", {
  p <- otpyr(a01())
  # a 2025: rows W/W, E/W, W/W; sorted E, W, W: within 1, across 1;
  # configurations 100, 000: ac 2; placements 6; cluster "East + West"
  expect_identical(as.character(p$id), c("a", "b", "a", "c"))
  expect_identical(p$vm, c(2L, 0L, 1L, 0L))
  expect_identical(p$cx, c(1L, 0L, 0L, 0L))
  expect_identical(p$np, c(6, 5, 4, 2))
  expect_identical(p$rc, c("East + West", "East", "East + North", "East"))
  expect_identical(attr(p, "years"), c(2025, 2026))
  # Ontario's own spelling of the columns is found
  d <- a01()
  names(d) <- c(
    "end_fiscal_year",
    "unique_individual_id",
    "region_at_time_of_placement",
    "region_most_recent_placement",
    "gender",
    "age_category",
    "mental_health_alert",
    "suicide_risk_alert",
    "suicide_watch_alert",
    "number_of_placements"
  )
  expect_identical(otpyr(d)$vm, p$vm)
  # the years kept; a missing value dropped and counted
  expect_identical(nrow(otpyr(a01(), years = 2026)), 2L)
  d <- a01()
  d$Gender[2] <- NA
  expect_identical(attr(otpyr(d), "drop"), 1L)
  errs(
    "otpyr",
    list(
      list("no column", a01()[-1]),
      list("must be named", a01(), cols = "x"),
      list("data frame or a table key", 1:3)
    )
  )
})

test_that("otpyr() reads a table of the service by its key", {
  dd <- odfix()
  gz <- file.path(dd, "db1", "a01.csv.gz")
  con <- gzfile(gz, "w")
  utils::write.csv(a01(), con, row.names = FALSE)
  close(con)
  odman(
    dd,
    sprintf(
      r"({"datasets": [{"db": "db1", "table": "a01", "key": "otis/a01",
    "rows": 7, "columns": [], "bytes_gz": %.0f, "sha256": "%s"}]})",
      file.size(gz),
      .odhash(gz, "test")
    )
  )
  odrun(dd, expect_identical(otpyr("otis/a01")$vm, otpyr(a01())$vm))
})

test_that("otpyr() takes other column names, refuses empty years, prints", {
  d <- a01()
  names(d)[5] <- "Sex"
  expect_identical(otpyr(d, cols = c(gender = "Sex")), otpyr(a01()))
  errs("otpyr", list(list("in the years asked", a01(), years = 1999)))
  d <- a01()
  d$Gender[2] <- NA
  p <- otpyr(d)
  # four person-years in 2025 and 2026, the row with no gender dropped
  prt(
    p,
    "^otpyr: 4 person-years, fiscal years 2025-2026 \\(2\\), 1 rows dropped"
  )
  expect_output(print(p), "year py  cx  vm")
  # nothing dropped, nothing said of it
  expect_false(any(grepl("dropped", capture.output(print(otpyr(a01()))),
                        fixed = TRUE)))
})

test_that("a flag that is neither yes nor no drops its record everywhere", {
  d <- data.frame(
    EndFiscalYear = 2025, UniqueIndividual_ID = c("a", "a", "b", "c", "c"),
    Region_AtTimeOfPlacement = c("E", "W", "E", "N", "N"),
    Region_MostRecentPlacement = "E", Gender = "Male", Age_Category = "x",
    MentalHealth_Alert = c("Yes", "No", "Unknown", "No", "No"),
    SuicideRisk_Alert = "No", SuicideWatch_Alert = "No",
    Number_Of_Placements = 1, stringsAsFactors = FALSE
  )
  p <- otpyr(d)
  # person b's only record has no usable flag: b is not a person-year, and
  # nobody's alert counts are another's
  expect_identical(p$id, c("a", "c"))
  expect_identical(attr(p, "drop"), 1L)
  expect_identical(p$ac, c(2L, 1L))
  expect_identical(p$cx, c(1L, 0L))
  expect_identical(p$anyal, c(1L, 0L))
  expect_identical(p$n, c(2L, 2L))
  rest <- otpyr(d[d$UniqueIndividual_ID != "b", ])
  expect_identical(p, rest, ignore_attr = "drop")
  d$MentalHealth_Alert <- "Unknown"
  errs("otpyr", list(list("no row has every column present", d)))
  # a number that is not 0 or 1 is no flag: it is not cut to one
  f <- alrtcx(data.frame(a = c(0.6, 1, 0, 1.9)), c(1, 1, 2, 2))
  expect_identical(f$n, c(1L, 1L))
  expect_identical(f$ac, c(1L, 1L))
})

test_that("every default address of a service is the one of odtab()", {
  ns <- asNamespace("crimlab")
  d <- unlist(lapply(ls(ns), function(f) {
    v <- get(f, ns)
    if (is.function(v)) {
      vapply(Filter(is.call, as.list(formals(v))), deparse1, "")
    }
  }))
  expect_identical(unique(grep("://", d, fixed = TRUE, value = TRUE)),
                   deparse1(formals(odtab)$base))
  expect_identical(formals(otpyr)$base, formals(odtab)$base)
})
