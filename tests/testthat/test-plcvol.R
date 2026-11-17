# Placement volatility (plcvol()) and alert complexity (alrtcx()) on small
# made-up records: every expected value is worked out by hand in the
# comments.

test_that("plcvol() counts changes within and between records", {
  # person 1: rows (W,W), (E,W), (W,W). Within: (E,W) = 1. Given order
  # W, E, W: two changes; sorted by region E, W, W: one.
  ra <- c("W", "E", "W", "N", "N")
  rb <- c("W", "W", "W", "N", "N")
  g <- plcvol(ra, rb, c(1, 1, 1, 2, 2), ord = "given")
  expect_identical(g$within, c(1L, 0L))
  expect_identical(g$across, c(2L, 0L))
  expect_identical(g$vm, c(3L, 0L))
  r <- plcvol(ra, rb, c(1, 1, 1, 2, 2))
  expect_identical(r$vm, c(2L, 0L))
  # sorted by region, across is the distinct regions at placement less one
  expect_identical(r$across, c(1L, 0L))
  # the bound 2R - 1 is reached by alternating records
  b <- plcvol(c("A", "B", "A"), c("B", "A", "B"), rep(1, 3), ord = "given")
  expect_identical(b$vm, 5L)
  # an injective relabelling changes nothing; merging regions can
  m <- c(W = "x", E = "y", N = "z")
  expect_identical(plcvol(m[ra], m[rb], c(1, 1, 1, 2, 2))$vm, r$vm)
  expect_identical(
    plcvol(rep("x", 5), rep("x", 5), c(1, 1, 1, 2, 2))$vm,
    c(0L, 0L)
  )
  # rows in another order: within unchanged; region order unchanged
  o <- c(3, 1, 2, 5, 4)
  expect_identical(plcvol(ra[o], rb[o], c(1, 1, 1, 2, 2)[o])$vm, r$vm)
  # a list id: person and year
  py <- plcvol(ra, rb, list(id = c(1, 1, 1, 1, 1), yr = c(1, 1, 1, 2, 2)))
  expect_identical(names(py)[1:2], c("id", "yr"))
  expect_identical(py$vm, c(2L, 0L))
  errs(
    "plcvol",
    list(
      list("one value per record", ra, rb[-1], 1:5),
      list("one label per record", ra, rb, 1:2),
      list("all present", NA, NA, 1),
      list("must be one of", ra, rb, 1:5, ord = "x")
    )
  )
})

test_that("alrtcx() counts distinct configurations, not alerts at once", {
  f <- data.frame(
    mh = c("Yes", "No", "No", "Yes", "Yes", NA),
    sr = c("No", "No", "Yes", "Yes", "Yes", "No"),
    sw = c("No", "No", "No", "Yes", "Yes", "No")
  )
  a <- alrtcx(f, c(1, 1, 2, 3, 3, 4))
  # 1: 100 and 000, two configurations, though never two alerts at once;
  # 2: 010 alone; 3: 111 twice, one configuration though three alerts;
  # 4: a missing flag, dropped
  expect_identical(a$id, c(1, 2, 3))
  expect_identical(a$ac, c(2L, 1L, 1L))
  expect_identical(a$cx, c(1L, 0L, 0L))
  expect_identical(a$anyal, c(1L, 1L, 1L))
  # keep counts only the configurations given
  expect_identical(
    alrtcx(f, c(1, 1, 2, 3, 3, 4), keep = "100")$ac,
    c(1L, 0L, 0L)
  )
  # logical and 0/1 flags give the same; order of rows does not matter
  l <- data.frame(x = c(TRUE, FALSE, TRUE), y = c(0, 0, 1))
  expect_identical(alrtcx(l, c(1, 1, 1))$ac, 3L)
  expect_identical(alrtcx(l[3:1, ], c(1, 1, 1))$ac, 3L)
  errs(
    "alrtcx",
    list(
      list("equal length", list(a = 1:2, b = 1), 1:2),
      list("codes of 2 digits", l, 1:3, keep = "1"),
      list("every flag present", data.frame(x = NA), 1)
    )
  )
})

test_that("plcvol() and alrtcx() name unnamed ids, take a matrix, print", {
  ra <- c("W", "E", "W", "N", "N")
  rb <- c("W", "W", "W", "N", "N")
  # an unnamed list of labels: its columns are named id1, id2, or id alone
  u <- plcvol(ra, rb, list(c(1, 1, 1, 2, 2), rep(2025, 5)))
  expect_identical(names(u)[1:2], c("id1", "id2"))
  expect_identical(u$vm, c(2L, 0L))
  expect_identical(names(plcvol(ra, rb, list(c(1, 1, 1, 2, 2))))[1], "id")
  # person 1: vm 2 (within 1, across 1); person 2: 0
  prt(
    plcvol(ra, rb, c(1, 1, 1, 2, 2)),
    paste0(
      "^plcvol: 2 person-years, records in region order; ",
      "mean volatility 1 \\(within 0.5, across 0.5\\)"
    )
  )
  prt(plcvol(ra, rb, c(1, 1, 1, 2, 2), ord = "given"), "in given order")
  # a matrix of flags, numeric or text, counts as its data frame does
  m <- cbind(
    mh = c(1, 0, 0, 1, 1),
    sr = c(0, 0, 1, 1, 1),
    sw = c(0, 0, 0, 1, 1)
  )
  id <- c(1, 1, 2, 3, 3)
  a <- alrtcx(m, id)
  expect_identical(a, alrtcx(as.data.frame(m), id))
  expect_identical(alrtcx(ifelse(m == 1, "Yes", "No"), id), a)
  expect_identical(a$ac, c(2L, 1L, 1L))
  # one of three person-years spans two configurations: 33.3%
  prt(
    a,
    paste0(
      "^alrtcx: 3 person-years; 1 \\(33.3%\\) span two or more ",
      "alert configurations"
    )
  )
})

test_that("labels holding a carriage return keep person-years apart", {
  id <- list(p = c("x\ry", "x"), y = c("z", "y\rz"))
  expect_identical(nrow(plcvol(c("A", "B"), c("A", "B"), id)), 2L)
  expect_identical(nrow(alrtcx(data.frame(a = c(1, 0)), id)), 2L)
})
