# Categorical integrity (R/catg.R): the tests of the same functions in
# rmoriebricklayer, carried over, then values computed once by
# rmoriebricklayer 1.1.1 on the same simulated input.

test_that("grdrcd maps by name, refuses unmapped, records the mapping", {
  r <- grdrcd(c("W", "B", "O", "W"), c(W = "White", B = "Black", O = "Other"))
  expect_identical(as.character(r), c("White", "Black", "Other", "White"))
  expect_error(grdrcd(c("W", "X"), c(W = "White")), "values with no mapping")
  expect_identical(
    as.character(grdrcd(c("W", "X", NA), c(W = "White"), keep = "X")),
    c("White", "X", NA)
  )
})

test_that("deccodes turns codes into labels and refuses unlabelled codes", {
  d <- deccodes(
    c(1, 2, 2, 3),
    c("1" = "White", "2" = "Black", "3" = "Indigenous")
  )
  expect_identical(as.character(d), c("White", "Black", "Black", "Indigenous"))
  expect_error(deccodes(c(1, 4), c("1" = "White")), "no mapping")
})

test_that("grdlvls declares levels and the reference", {
  f <- grdlvls(c("White", "Black", "White"), c("White", "Black"), "White")
  expect_identical(as.integer(f), c(1L, 2L, 1L))
  expect_error(
    grdlvls(c("White", "Other"), c("White", "Black")),
    "outside the declared levels"
  )
  expect_error(grdlvls("White", c("White", "Black"), "Black"),
               "is not the first level")
})


test_that("vrfrcd proves a recode and names a swap", {
  expect_silent(vrfrcd(
    c("W", "B", "W"),
    c("White", "Black", "White"),
    c(W = "White", B = "Black")
  ))
  expect_error(
    vrfrcd(c("W", "B"), c("Black", "White"), c(W = "White", B = "Black")),
    "THIS is how groups get swapped"
  )
  expect_error(vrfrcd(c("W", "B"), "White", c(W = "White")), "rows were lost")
})

test_that("vrfmrg accepts matching counts and names a permutation", {
  x <- c("White", "White", "Black", "Indigenous", "White", "Black")
  expect_true(vrfmrg(x, c(White = 3, Black = 2, Indigenous = 1))$ok)
  expect_error(
    vrfmrg(x, c(White = 2, Black = 3, Indigenous = 1)),
    "White -> Black"
  )
})


test_that("grdbin refuses categorical treatment columns", {
  expect_true(grdbin(c(0, 1, 1, 0), "d"))
  expect_error(grdbin(c("treated", "control"), "d"), "level indices")
  expect_error(grdbin(c(0, 2), "d"), "must be 0 or 1")
})


test_that("relabel refuses positional labels and maps by name", {
  f <- factor(c("W", "B", "O", "W"))
  expect_error(relbl(f, c("Black", "Other", "White")), "position")
  g <- relbl(f, c(W = "White", B = "Black", O = "Other"))
  expect_identical(as.character(g), c("White", "Black", "Other", "White"))
  expect_identical(levels(g), c("White", "Black", "Other"))
  expect_error(relbl(f, c(W = "White")), "no mapping")
})

test_that("declbld decodes by code, levels in code order", {
  x <- structure(
    c(1, 2, 2, 4),
    labels = c(White = 1, Black = 2, Other = 3, Unknown = 4)
  )
  f <- declbld(x)
  expect_identical(as.character(f), c("White", "Black", "Black", "Unknown"))
  expect_identical(levels(f), c("White", "Black", "Other", "Unknown"))
  expect_error(declbld(c(1, 5), c("1" = "White")), "no mapping")
  expect_error(declbld(c(1, 2)), "not the categories")
})


test_that("xfrvrf accepts a faithful import, names a rotated one", {
  x <- structure(
    c(1, 1, 2, 4, 1),
    labels = c(White = 1, Black = 2, Other = 3, Unknown = 4)
  )
  ok <- xfrvrf(
    x,
    c(White = 3, Black = 1, Unknown = 1),
    codebook = c("1" = "White", "2" = "Black", "3" = "Other", "4" = "Unknown")
  )
  expect_true(ok$ok)
  expect_true(ok$codebook)
  expect_identical(levels(ok$decoded), c("White", "Black", "Other", "Unknown"))
  rotated <- structure(
    c(1, 1, 2, 4, 1),
    labels = c(Black = 1, Other = 2, Unknown = 3, White = 4)
  )
  expect_error(
    xfrvrf(
      rotated,
      c(White = 3, Black = 1, Unknown = 1),
      codebook = c("1" = "White", "2" = "Black", "3" = "Other", "4" = "Unknown")
    ),
    "disagree with the source code book"
  )
  expect_error(
    xfrvrf(
      c("Black", "Black", "Other", "White", "Black"),
      c(White = 3, Black = 1, Unknown = 1)
    ),
    "do not match|permuted|not in the published"
  )
})

test_that("deccodes keeps or drops NA codes as asked", {
  d <- deccodes(c(1, NA, 2), c("1" = "White", "2" = "Black"))
  expect_identical(as.character(d), c("White", NA, "Black"))
  expect_error(deccodes(c(1, 2), c("White", "Black")), "named character")
})

test_that("grdrcd refuses an unnamed mapping and factor input works", {
  expect_error(
    grdrcd(c("W", "B"), c("White", "Black")),
    "named character vector"
  )
  r <- grdrcd(factor(c("W", "B")), c(W = "White", B = "Black"))
  expect_identical(as.character(r), c("White", "Black"))
})


test_that("vrfrcd reports fan-out and missingness leaks", {
  expect_error(
    vrfrcd(c("W", "W"), c("White", "Black"), c(W = "White")),
    "MULTIPLE|fan"
  )
  expect_error(
    vrfrcd(c("W", "B"), c("White", NA), c(W = "White", B = "Black")),
    "missing"
  )
})

test_that("vrfmrg honours a tolerance, reports mismatches", {
  x <- c("White", "White", "Black")
  expect_true(vrfmrg(x, c(White = 3, Black = 1), tol = 1)$ok)
  expect_error(vrfmrg(x, c(White = 5, Black = 9)), "do not match")
  expect_error(vrfmrg(x, c(3, 1)), "named numeric")
  # nothing published compares nothing: all() of no categories is not a pass
  empty <- stats::setNames(numeric(), character())
  expect_error(vrfmrg(rep(NA, 3), empty, strict = FALSE), "named numeric")
  expect_error(vrfmrg(x, c(White = NA, Black = 1)), "named numeric")
})


test_that("grdbin accepts NA, refuses factors and non-binary numbers", {
  expect_true(grdbin(c(0, 1, NA), "d"))
  expect_error(grdbin(factor(c("a", "b")), "d"), "level indices")
  expect_error(grdbin(c(0.5, 1), "d"), "0 or 1")
})


test_that("relabel keeps pass-through levels and refuses unmapped ones", {
  f <- factor(c("W", "B", "X"))
  g <- relbl(f, c(W = "White", B = "Black"), keep = "X")
  expect_identical(as.character(g), c("White", "Black", "X"))
  expect_true("X" %in% levels(g))
  expect_error(relbl(c("W", "Q"), c(W = "White")), "no mapping")
})

test_that("declbld takes explicit labels, refuses bad attributes", {
  d <- declbld(c(2, 1), vlabs = c("1" = "White", "2" = "Black"))
  expect_identical(as.character(d), c("Black", "White"))
  expect_identical(levels(d), c("White", "Black"))
  bad <- structure(c(1, 2), labels = c(1, 2))
  expect_error(declbld(bad), "no names")
  s <- declbld(c("b", "a"), vlabs = c(a = "Alpha", b = "Beta"))
  expect_identical(levels(s), c("Alpha", "Beta"))
})


test_that("xfrvrf decodes plain codes with a code book, tolerance", {
  cb <- c("1" = "White", "2" = "Black")
  r <- xfrvrf(c(1, 1, 2), c(White = 2, Black = 1), codebook = cb)
  expect_true(r$ok)
  expect_true(is.na(r$codebook))
  expect_identical(levels(r$decoded), c("White", "Black"))
  lab <- structure(c(1, 1, 2), labels = c(White = 1, Black = 2))
  t2 <- xfrvrf(lab, c(White = 3, Black = 1), tol = 1)
  expect_true(t2$ok)
  expect_error(
    xfrvrf(c("White", "Other"), c(White = 1, Black = 1)),
    "not in the published"
  )
})


test_that("vrfmrg and xfrvrf return ok = FALSE when asked", {
  x <- c("White", "White", "Black", "Indigenous", "White", "Black")
  r <- vrfmrg(x, c(White = 2, Black = 3, Indigenous = 1), strict = FALSE)
  expect_false(r$ok)
  expect_identical(unname(r$permutation[["White"]]), "Black")
  expect_match(r$message, "permuted")
  r2 <- vrfmrg(
    c(x, "Other"),
    c(White = 3, Black = 2, Indigenous = 1),
    strict = FALSE
  )
  expect_false(r2$ok)
  expect_match(r2$message, "not in the published")
  cb <- c("1" = "White", "2" = "Black", "3" = "Other", "4" = "Unknown")
  rotated <- structure(
    c(1, 1, 2, 4, 1),
    labels = c(Black = 1, Other = 2, Unknown = 3, White = 4)
  )
  t <- xfrvrf(
    rotated,
    c(White = 3, Black = 1, Unknown = 1),
    codebook = cb,
    strict = FALSE
  )
  expect_false(t$ok)
  expect_false(t$codebook)
  expect_true(any(grepl("code book", t$reasons, fixed = TRUE)))
  # the forensics input is reachable from the public path now
  t2 <- xfrvrf(
    c("Black", "White", "White", "Black", "Black"),
    c(White = 3, Black = 2),
    strict = FALSE
  )
  expect_false(t2$ok)
  expect_identical(unname(t2$marginals$permutation[["White"]]), "Black")
  expect_error(
    xfrvrf(
      c("Black", "White", "White", "Black", "Black"),
      c(White = 3, Black = 2)
    ),
    "do not match"
  )
})



test_that("the results equal rmoriebricklayer's", {
  # rmoriebricklayer 1.1.1 on the same simulated input, printed with dput()
  m <- c("1" = "White", "2" = "Black", "3" = "Indigenous")
  expect_identical(as.vector(deccodes(c(3, 1, NA, 2, 2), m)),
                   c("Indigenous", "White", NA, "Black", "Black"))
  r <- relbl(c("b", "a", "c", "a"), c(c = "C", a = "A", b = "B"))
  expect_identical(levels(r), c("C", "A", "B"))
  expect_identical(as.integer(r), c(3L, 2L, 1L, 2L))
  r <- relbl(factor(c("x", "y"), levels = c("y", "x", "z")),
             c(y = "Y", x = "X"), keep = "z")
  expect_identical(levels(r), c("Y", "X"))
  expect_identical(as.integer(r), 2:1)
  d <- declbld(c(10, 2, 2, 10), c("10" = "Ten", "2" = "Two"))
  expect_identical(levels(d), c("Two", "Ten"))
  expect_identical(as.integer(d), c(2L, 1L, 1L, 2L))
  expect_identical(grdlvls(c("b", "a"), c("b", "a"), "b"),
                   factor(c("b", "a"), levels = c("b", "a")))
  v <- vrfrcd(c("M", "F", "F", NA, "X"), c("Male", "Female", "Female", NA, "X"),
              c(M = "Male", F = "Female"))
  expect_identical(as.character(v$original), rep(c("F", "M", "X"), 3))
  expect_identical(as.character(v$recoded),
                   rep(c("Female", "Male", "X"), each = 3))
  expect_identical(v$Freq, c(2L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L))
  x <- rep(c("White", "Black", "Indigenous"), c(50, 30, 20))
  r <- vrfmrg(x, c(White = 30, Black = 50, Indigenous = 20), strict = FALSE)
  expect_identical(r$counts, c(White = 50, Black = 30, Indigenous = 20))
  expect_false(r$ok)
  expect_identical(r$permutation, c(White = "Black", Black = "White",
                                    Indigenous = "Indigenous"))
  expect_identical(r$message, paste0(
    "vrfmrg: recoded counts do not match the published counts. White: ",
    "observed 50, published 30; Black: observed 30, published 50; ",
    "Indigenous: observed 20, published 20. The observed counts match the ",
    "published ones if the labels are permuted (White -> Black, Black -> ",
    "White): the labels are attached to the wrong groups."))
  expect_output(print(r), "the counts differ")
  r <- vrfmrg(x, c(White = 50, Black = 30, Indigenous = 21), tol = 1)
  expect_true(r$ok)
  expect_null(r$permutation)
  y <- xfrvrf(c(1, 2, 2, 1, 1), c(No = 3, Yes = 2),
              codebook = c("1" = "No", "2" = "Yes"))
  expect_true(y$ok)
  expect_identical(as.integer(y$decoded), c(1L, 2L, 2L, 1L, 1L))
  expect_identical(levels(y$decoded), c("No", "Yes"))
  expect_output(print(y), "as the source reported")
})
