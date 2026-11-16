# Region-map checks (R/rgnmap.R): the tests of the same functions in
# rmoriebricklayer 1.1.1, carried over; the expected values there were
# written for rmoriebricklayer and are the values it gives.

test_that("rgnmapintg finds duplicates, gaps and unknown regions", {
  cw <- data.frame(
    inst = c("North", "South", "East"),
    cd = c("3557", "3520", "3506"),
    stringsAsFactors = FALSE
  )
  known <- c("3557", "3520", "3506")
  ok <- rgnmapintg(cw, "inst", "cd", regs = known)
  expect_true(all(ok$pass))
  expect_identical(nrow(ok), 3L)

  dup <- rbind(cw, cw[1, ])
  expect_identical(rgnmapintg(dup, "inst", "cd", known)$observed[1], 1)

  gap <- cw
  gap$cd[2] <- ""
  expect_identical(rgnmapintg(gap, "inst", "cd", known)$observed[2], 1)

  wrong <- cw
  wrong$cd[3] <- "2406"
  expect_identical(rgnmapintg(wrong, "inst", "cd", known)$observed[3], 1)

  ## without a reference geography the third check is not offered at all,
  ## rather than passing vacuously
  expect_identical(nrow(rgnmapintg(cw, "inst", "cd")), 2L)
})

test_that("rgnmapintg errors on a column that is not there", {
  cw <- data.frame(inst = "North", cd = "3557", stringsAsFactors = FALSE)
  expect_error(rgnmapintg(cw, "facility", "cd"), "no column")
  expect_error(rgnmapintg(cw, "inst", "division"), "no column")
})

test_that("rgnmapcmp reports agreement and the first disagreement", {
  pub <- data.frame(
    inst = c("North", "South"),
    cd = c("3557", "3520"),
    stringsAsFactors = FALSE
  )
  same <- rgnmapcmp(pub, pub, "inst")
  expect_true(all(same$mismatched == 0))
  expect_true(all(same$first == ""))

  obs <- pub
  obs$cd[2] <- "3521"
  diff <- rgnmapcmp(pub, obs, "inst")
  cd <- diff[diff$column == "cd", ]
  expect_equal(cd$mismatched, 1, tolerance = testthat::testthat_tolerance())
  expect_match(cd$first, "South")
  expect_match(cd$first, "3521")
})

test_that("rgnmapcmp counts units present on one side only", {
  pub <- data.frame(
    inst = c("North", "South"),
    cd = c("3557", "3520"),
    stringsAsFactors = FALSE
  )
  obs <- pub[1, ]
  rows <- rgnmapcmp(pub, obs, "inst")
  rows <- rows[rows$column == "rows", ]
  expect_equal(rows$mismatched, 1, tolerance = testthat::testthat_tolerance())
  expect_match(rows$first, "published only")
})

test_that("rgnmapcmp tolerates a numeric round trip through text", {
  pub <- data.frame(inst = "North", lat = 46.55004065, stringsAsFactors = FALSE)
  obs <- data.frame(
    inst = "North",
    lat = as.numeric(format(46.55004065, digits = 15)),
    stringsAsFactors = FALSE
  )
  cmp <- rgnmapcmp(pub, obs, "inst", cols = "lat")
  expect_equal(
    cmp$mismatched[cmp$column == "lat"],
    0,
    tolerance = testthat::testthat_tolerance()
  )
})

test_that("rgnmapalt separates documented from unexplained", {
  cw <- data.frame(
    inst = c("North", "South", "Hill"),
    cd = c("3557", "3520", "3506"),
    stringsAsFactors = FALSE
  )
  route <- c("North" = "3557", "South" = "3520", "Hill" = "3519")

  d <- rgnmapalt(cw, "inst", "cd", route, known = "Hill")
  expect_identical(nrow(d), 1L)
  expect_true(d$known)
  expect_identical(sum(!d$known), 0L)

  route["South"] <- "3521"
  d <- rgnmapalt(cw, "inst", "cd", route, known = "Hill")
  expect_identical(nrow(d), 2L)
  expect_identical(sum(!d$known), 1L)
  expect_identical(d$unit[!d$known], "South")
})

test_that("rgnmapalt skips units the route does not cover", {
  cw <- data.frame(
    inst = c("North", "South"),
    cd = c("3557", "3520"),
    stringsAsFactors = FALSE
  )
  ## the route knows nothing about South, which is not a disagreement
  d <- rgnmapalt(cw, "inst", "cd", c("North" = "3557"))
  expect_identical(nrow(d), 0L)
})

test_that("rgnmapalt accepts the route as a data frame", {
  cw <- data.frame(
    inst = c("North", "South"),
    cd = c("3557", "3520"),
    stringsAsFactors = FALSE
  )
  rt <- data.frame(
    inst = c("North", "South"),
    cd = c("3557", "3599"),
    stringsAsFactors = FALSE
  )
  d <- rgnmapalt(cw, "inst", "cd", rt)
  expect_identical(d$unit, "South")
  expect_identical(d$second, "3599")
})


test_that("the results print and refuse what they cannot read", {
  cw <- data.frame(inst = c("A", "B"), cd = c("1", "2"),
                   stringsAsFactors = FALSE)
  expect_s3_class(rgnmapintg(cw, "inst", "cd"), "rgnmapintg")
  expect_output(print(rgnmapintg(cw, "inst", "cd")), "every check passed")
  expect_output(print(rgnmapintg(rbind(cw, cw), "inst", "cd")), "1 of 2")
  expect_output(print(rgnmapcmp(cw, cw, "inst")), "the maps agree")
  o <- cw
  o$cd[1] <- "9"
  expect_output(print(rgnmapcmp(cw, o, "inst")), "the maps disagree")
  expect_error(rgnmapcmp(cw, o, "inst", cols = "x"), "missing: x")
  expect_error(rgnmapintg(list(), "inst", "cd"), "must be a data frame")
  expect_error(rgnmapalt(cw, "inst", "cd", c("1", "2")), "named by unit")
  expect_identical(nrow(rgnmapalt(cw, "inst", "cd", c(A = "1", B = NA))), 0L)
})

test_that("rgnmapcmp names a unit found only in the recomputed map", {
  pub <- data.frame(inst = "North", cd = "1", stringsAsFactors = FALSE)
  obs <- data.frame(inst = c("North", "South"), cd = c("1", "2"),
                    stringsAsFactors = FALSE)
  r <- rgnmapcmp(pub, obs, "inst")
  expect_identical(r$first[r$column == "rows"], "recomputed only: South")
  expect_identical(r$mismatched[r$column == "rows"], 1L)
})

test_that("arguments that are not names are refused in words", {
  m <- data.frame(inst = c("a", "b"), cd = c("1", "2"))
  bad <- rawToChar(as.raw(c(0xff, 0xfe, 0x41)))
  expect_error(expect_warning(rgnmapalt(m, "inst", "cd", sum), NA),
               "^rgnmapalt: `route` must be a vector of names or codes")
  expect_error(rgnmapalt(m, "inst", "cd", c(a = "1"), known = sum),
               "^rgnmapalt: `known` must be a vector of names or codes")
  expect_error(rgnmapalt(m, sum, "cd", c(a = "1")),
               "^rgnmapalt: `unit` must be one column name of `map`")
  expect_error(rgnmapalt(m, "inst", "cd", c(a = bad)),
               "^rgnmapalt: `route` holds text that is not valid")
  expect_error(rgnmapintg(m, "inst", "cd", regs = new.env()),
               "^rgnmapintg: `regs` must be a vector")
  expect_error(rgnmapcmp(m, m, "inst", cols = sum),
               "^rgnmapcmp: `cols` must be a vector")
})
