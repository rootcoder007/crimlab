# catchk(): counts against table(), each problem it reports, the reference
# it names and why, the stated reference, the input forms and the refusals.

lv <- function(k, v) k$levels[k$levels$var == v, ]
vr <- function(k, v) k$vars[k$vars$var == v, ]

test_that("counts, missing values and levels agree with table()", {
  f <- factor(c("b", NA, "a", "b", "b"), levels = c("a", "b", "c"))
  k <- catchk(f, nmin = 2)
  l <- lv(k, "x")
  expect_identical(l$level, levels(f))
  expect_identical(l$n, as.vector(table(f), "integer"))
  v <- vr(k, "x")
  expect_identical(v$type, "factor")
  expect_identical(v$n, 5L)
  expect_identical(v$nna, sum(is.na(f)))
  expect_identical(v$nlev, 3L)
  expect_identical(v$unused, 1L)
  expect_identical(v$low, 1L)
  expect_identical(l$low, c(TRUE, FALSE, FALSE))
  expect_identical(v$ref, "a")
  expect_identical(v$why, "factor: the first level")
  expect_identical(l$ref, c(TRUE, FALSE, FALSE))
  expect_false(v$ok)
  expect_false(k$ok)
  # character: first appearance, no default reference
  s <- c("z", "y", NA, "z", "x", "y", "z")
  k <- catchk(s, nmin = 1)
  l <- lv(k, "x")
  expect_identical(l$level, c("z", "y", "x"))
  expect_identical(l$n, as.vector(table(s)[l$level], "integer"))
  expect_identical(vr(k, "x")$ref, NA_character_)
  expect_identical(vr(k, "x")$why, "no order: state it")
  expect_identical(vr(k, "x")$nna, 1L)
  expect_true(catchk(c("a", "b", "a"), nmin = 1)$ok)
  expect_identical(catchk(ordered(c("lo", "hi")), nmin = 1)$vars$type,
                   "ordered")
})

test_that("logical and 0/1 codings list both values and name 0 or FALSE", {
  k <- catchk(list(l = c(TRUE, FALSE, NA, TRUE), b = c(1, 0, 1, 1),
                   i = c(0L, 1L, NA), t = c(TRUE, TRUE)), nmin = 1)
  expect_identical(k$vars$bin, c(TRUE, TRUE, TRUE, TRUE))
  expect_identical(k$vars$ref, c("FALSE", "0", "0", "FALSE"))
  expect_identical(k$vars$why[1:2], c("logical: TRUE is the treatment",
                                      "0 and 1: 1 is the treatment"))
  expect_identical(lv(k, "b")$n, c(1L, 3L))
  expect_identical(lv(k, "l")$level, c("FALSE", "TRUE"))
  expect_identical(lv(k, "i")$n, c(1L, 1L))
  # only TRUE seen: FALSE is unused and fewer than two levels are seen
  expect_identical(vr(k, "t")$unused, 1L)
  expect_identical(k$vars$ok, c(FALSE, TRUE, FALSE, FALSE))
  expect_output(print(k),
                "t: never seen: \"FALSE\"; fewer than two levels seen")
})

test_that("labels that differ only by case or white space are twins", {
  s <- c("Male", " male", "MALE ", "Female", "female  adult",
         "Female Adult", "Été", "été")
  k <- catchk(s, nmin = 1)
  l <- lv(k, "x")
  expect_identical(vr(k, "x")$twin, 5L)
  expect_identical(l$twin, c(" male", "Male", "Male", NA, "Female Adult",
                             "female  adult", NA, NA))
  expect_false(k$ok)
  expect_output(print(k), paste0("differ only by case or spaces: \"Male\", ",
                                 "\" male\", \"MALE \"; \"female  adult\", ",
                                 "\"Female Adult\""), fixed = TRUE)
  # white space is made one space, never removed: "a b" and "ab" differ
  expect_identical(catchk(c("a b", "ab"), nmin = 1)$vars$twin, 0L)
})

test_that("whole numbers with few values look like category codes", {
  k <- catchk(list(a = c(1L, 2L, 2L, 1L), b = c(1, 5, 7), c = 1:11,
                   d = c(0.5, 1), e = c(1, Inf), f = c(2, NaN, 2)),
              nmin = 1)
  expect_identical(k$vars$codes, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(k$vars$type, c("integer", "double", "integer", "double",
                                  "double", "double"))
  expect_identical(k$vars$nlev, c(2L, 3L, 11L, NA, NA, 1L))
  expect_identical(k$vars$why[4:5], c("not categorical", "not categorical"))
  expect_identical(k$vars$nna, c(0L, 0L, 0L, 0L, 0L, 1L))
  expect_identical(k$vars$ok, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(nrow(lv(k, "d")), 0L)
})

test_that("value labels are checked against the values", {
  lb <- structure(c(1, 2, 9, 1, NA),
                  labels = setNames(1:3 + 0, c("Male", "Female", "Male")),
                  class = c("haven_labelled", "vctrs_vctr", "double"))
  k <- catchk(lb, nmin = 1)
  v <- k$vars
  expect_identical(v$type, "labelled")
  expect_identical(v$nolab, 1L)
  expect_identical(v$dlab, 2L)
  expect_identical(v$unused, 1L)
  expect_identical(v$nna, 1L)
  expect_identical(v$ref, NA_character_)
  expect_identical(v$why, "labelled: decode it first")
  expect_false(v$bin)
  expect_identical(k$levels$level, c("1", "2", "9", "3"))
  expect_identical(k$levels$label, c("Male", "Female", NA, "Male"))
  expect_identical(k$levels$n, c(2L, 1L, 1L, 0L))
  expect_output(print(k), "value labels: decode it first")
  # labels with no names, and a labelled character and logical vector
  k <- catchk(list(a = structure(c(1, 2), labels = c(1, 2)),
                   b = structure(c("m", "f"), labels = c(Man = "m")),
                   c = structure(c(TRUE, FALSE), labels = c(yes = TRUE))),
              nmin = 1)
  expect_identical(k$vars$nolab, c(2L, 1L, 1L))
  expect_identical(k$vars$bin, c(FALSE, FALSE, FALSE))
  expect_identical(lv(k, "b")$label, c("Man", NA))
  expect_identical(k$vars$ok, c(FALSE, FALSE, FALSE))
})

test_that("variables that are not categorical are reported, not audited", {
  d <- data.frame(day = as.Date("2026-01-01") + 0:1,
                  t = as.POSIXct("2026-01-01", tz = "UTC") + 0:1,
                  w = c(0.25, NA))
  d$dt <- as.difftime(c(1, 2), units = "days")
  d$m <- matrix(1:4, 2)
  d$l <- list(1, "a")
  k <- catchk(d)
  expect_identical(k$vars$type, c("Date", "POSIXct", "double", "difftime",
                                  "matrix", "list"))
  expect_identical(k$vars$nna, c(0L, 0L, 1L, 0L, 0L, 0L))
  expect_identical(k$vars$ok, c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(nrow(k$levels), 0L)
  expect_identical(k$vars$state, rep(NA_character_, 6))
  k <- catchk(list(lm = matrix(c(TRUE, NA), 1), cm = matrix(c("a", NA), 1),
                   z = complex(real = 1:2), r = as.raw(1:2)))
  expect_identical(k$vars$type, c("matrix", "matrix", "complex", "raw"))
  expect_identical(k$vars$nna, c(1L, 1L, 0L, 0L))
})

test_that("a stated reference is checked against the levels", {
  f <- factor(c("b", "a", "b", "B", "b", "a", "a"), levels = c("a", "b", "B",
                                                               "c"))
  st <- function(r, nmin = 2) catchk(f, ref = r, nmin = nmin)$vars$state
  expect_identical(st("a"), "default")
  expect_identical(st("b"), "other")
  expect_identical(st("B"), "low")
  expect_identical(st("c"), "empty")
  expect_identical(st("A"), "case")
  expect_identical(st(" b "), "case")
  expect_identical(st("d"), "absent")
  expect_identical(st(factor("b")), "other")
  # by variable name, a list, and a number read as text
  k <- catchk(list(f = f, n = c(1, 2, 2, 1)), ref = list(n = 2, f = "a"),
              nmin = 1)
  expect_identical(k$vars$want, c("a", "2"))
  expect_identical(k$vars$state, c("default", "other"))
  k <- catchk(list(f = f, n = c(1, 2)), ref = c(n = "3"), nmin = 1)
  expect_identical(k$vars$state, c(NA, "absent"))
  expect_false(k$vars$ok[2])
  # one level for every variable; a variable that is not categorical
  k <- catchk(data.frame(a = c("x", "y"), b = c(0.5, 1)), ref = "x",
              nmin = 1)
  expect_identical(k$vars$state, c("other", "absent"))
  expect_output(print(catchk(f, ref = "d")), "reference \"d\": not a level")
  expect_output(print(catchk(f, ref = "c")), "reference \"c\": never seen")
  expect_output(print(catchk(f, ref = "A")),
                "reference \"A\": present only with other case or spacing")
  expect_output(print(catchk(f, ref = "B", nmin = 2)),
                "reference \"B\": seen fewer than nmin times")
})

test_that("the input forms name the variables", {
  expect_identical(catchk(list(1, 2))$vars$var, c("1", "2"))
  expect_identical(catchk(list(a = 1, 2))$vars$var, c("a", "2"))
  expect_identical(catchk(c(a = 1, b = 2))$vars$var, "x")
  expect_identical(catchk(data.frame(p = 1, q = "a"))$vars$var, c("p", "q"))
  e <- catchk(character(0))
  expect_identical(e$vars$nlev, 0L)
  expect_false(e$ok)
  expect_identical(catchk(factor(character(0), levels = character(0)))$vars$why,
                   "factor with no levels")
  expect_identical(catchk(list(NULL))$vars$type, "NULL")
})

test_that("catchk() refuses what it cannot read", {
  f <- c("a", "b")
  errs("catchk", list(
    list("must be a vector, a data frame or a list", NULL),
    list("must be a vector, a data frame or a list", mean),
    list("not a matrix", matrix(1:4, 2)),
    list("has no variables", list()),
    list("names a variable twice: a", setNames(list(1, 2), c("a", "a"))),
    list("`nmin` must be one whole number from 1", f, NULL, 0),
    list("`nmin` must be one whole number from 1", f, NULL, 1.5),
    list("`ref` must be levels", f, mean),
    list("`ref` must be levels", f, Sys.Date()),
    list("`ref` must name at least one level", f, character(0)),
    list("with no missing values", f, NA_character_),
    list("with no missing values", f, list(1:2)),
    list("one level, or levels named by variable", f, c("a", "b")),
    list("`ref` names no variable: y", f, c(y = "a")),
    list("`ref` names a variable twice: x", f,
         setNames(c("a", "b"), c("x", "x"))),
    list("has value labels", f, structure(1, labels = c(a = 1)))
  ))
})

test_that("print.catchk() lists the problems of each variable", {
  d <- data.frame(sex = c("Male", "male", "Female", "Female", NA),
                  custody = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                  court = c(1L, 2L, 2L, 3L, 1L))
  k <- catchk(d, nmin = 2)
  out <- capture.output(print(k))
  expect_identical(out[1], "catchk: 3 variables, 2 with problems")
  expect_true(any(out == paste0("sex: 1 missing; seen fewer than 2 times: ",
                                "\"Male\", \"male\"; differ only by case or ",
                                "spaces: \"Male\", \"male\"")))
  expect_true(any(out == paste0("court: seen fewer than 2 times: \"3\"; ",
                                "numbers that look like category codes")))
  expect_output(print(catchk(c(1, 2), nmin = 1)), "catchk: 1 variable, 1 with")
  two <- structure(c(1, 2), labels = setNames(c(1, 2), c("a", "a")))
  expect_output(print(catchk(two, nmin = 1)), "2 labels shared by two values")
  ok <- capture.output(print(catchk(c("a", "b"), nmin = 1)))
  expect_length(ok, 3L)
  expect_output(expect_invisible(print(k)))
})
