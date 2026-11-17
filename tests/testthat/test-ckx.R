# The argument checks of R/ckx.R, shared by the search and rescue, air
# pollution, footprint and data description families: each refusal names
# the function and the argument.

test_that("a refusal names the function and ends in a full stop", {
  expect_error(.ckstop("f", "`%s` is wrong", "x"), "^f: `x` is wrong\\.$")
})

test_that("a number is kept as a double and refused for each rule", {
  expect_identical(.cknum(2L, "n", "f"), 2)
  expect_identical(dim(.cknum(matrix(1:4, 2), "m", "f", scalar = FALSE)),
                   c(2L, 2L))
  expect_error(.cknum(NULL, "n", "f"), "f: `n` must be numeric, got NULL.",
               fixed = TRUE)
  expect_error(.cknum(TRUE, "n", "f"), "not logical")
  expect_error(.cknum("1", "n", "f"), "got character")
  expect_error(.cknum(numeric(), "n", "f"), "at least one value")
  expect_error(.cknum(1:2, "n", "f"), "a single number")
  expect_error(.cknum(Inf, "n", "f"), "finite, with no missing values")
  expect_error(.cknum(0, "n", "f", positive = TRUE), "must be positive")
  expect_error(.cknum(5, "n", "f", min = 0, max = 1), "between 0 and 1")
  expect_error(.cknum(-1, "n", "f", min = 0), "at least 0")
  expect_error(.cknum(1.5, "n", "f", whole = TRUE), "a whole number")
})

test_that("a choice is matched partially, the first by default", {
  ch <- c("point", "line")
  expect_identical(.ckarg(NULL, ch, "w", "f"), "point")
  expect_identical(.ckarg(ch, ch, "w", "f"), "point")
  expect_identical(.ckarg("li", ch, "w", "f"), "line")
  expect_error(.ckarg("z", ch, "w", "f"),
               "f: `w` must be one of \"point\", \"line\", not \"z\".",
               fixed = TRUE)
})

test_that("flags, text and seeds are checked", {
  expect_true(.ckflag(TRUE, "q", "f"))
  expect_error(.ckflag(NA, "q", "f"), "`q` must be TRUE or FALSE")
  expect_identical(.cktxt("abc", "s", "f"), "abc")
  expect_identical(.ckseed(3, "seed", "f"), 3L)
  expect_identical(.ckseed(2^40, "seed", "f"), 2^40)
  expect_error(.ckseed(1.5, "seed", "f"), "single whole number")
})

test_that("a seeded expression leaves the caller's stream as it was", {
  set.seed(1)
  s <- get(".Random.seed", envir = globalenv())
  a <- .clseed(5L, stats::runif(2))
  expect_identical(get(".Random.seed", envir = globalenv()), s)
  expect_identical(.clseed(5L, stats::runif(2)), a)
  b <- .clseed(2^40, stats::runif(2))
  expect_identical(.clseed(2^40, stats::runif(2)), b)
  expect_false(identical(.clseed(2^40 + 1, stats::runif(2)), b))
  expect_identical(get(".Random.seed", envir = globalenv()), s)
})

test_that(".clas() gives another function's refusal under the caller's name", {
  expect_error(.clas("runfx", stop("ppmlfe: `x` is bad.", call. = FALSE)),
               "^runfx: `x` is bad\\.$")
  expect_error(.clas("runfx", stop("runfx: own.", call. = FALSE)),
               "^runfx: own\\.$")
  expect_error(.clas("runfx", stop("nofn: kept.", call. = FALSE)),
               "^nofn: kept\\.$")
  expect_warning(.clas("runfx", warning("gplume: w.", call. = FALSE)),
                 "^runfx: w\\.$")
  expect_identical(.clas("runfx", 3), 3)
})

test_that(".clas() renames a message that quotes invalid text", {
  m <- paste0("declbld: values with no mapping: '",
              rawToChar(as.raw(c(0xff, 0x41))), "'.")
  e <- tryCatch(.clas("xfrvrf", stop(m, call. = FALSE)),
                error = conditionMessage)
  expect_identical(rawToChar(charToRaw(e)[1:8]), "xfrvrf: ")
})

test_that("a condition quoting invalid text carries it as <xx>", {
  bad <- rawToChar(as.raw(c(0xff, 0xfe)))
  ok <- function(e) validEnc(conditionMessage(e))
  # refusals made in R and in the C++ binding
  expect_identical(tryCatch(clprep(preset = bad), error = ok), TRUE)
  expect_identical(tryCatch(dastd(c(5, 3), c(100, 100), c("a", "b"),
                                  ref = bad), error = ok), TRUE)
  # the package's stop(), warning() and message(), and conditions given
  expect_identical(tryCatch(stop("x ", bad), error = ok), TRUE)
  expect_identical(tryCatch(stop(simpleError(bad)), error = ok), TRUE)
  expect_identical(tryCatch(warning("x ", bad), warning = ok), TRUE)
  expect_identical(tryCatch(warning(simpleWarning(bad)), warning = ok), TRUE)
  expect_identical(tryCatch(message("x ", bad), message = ok), TRUE)
  expect_identical(tryCatch(message(simpleMessage(bad)), message = ok), TRUE)
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]))
  expect_identical(tryCatch(stop("x ", bad), error = conditionMessage),
                   "x <ff><fe>")
})

test_that("the package's stop() and warning() keep the call as base R's", {
  f <- function() stop("a")
  expect_identical(conditionCall(tryCatch(f(), error = identity)), quote(f()))
  g <- function() stop("a", call. = FALSE)
  expect_null(conditionCall(tryCatch(g(), error = identity)))
  h <- function() warning("w")
  expect_identical(conditionCall(tryCatch(h(), warning = identity)),
                   quote(h()))
  expect_identical(tryCatch(message("m"), message = conditionMessage), "m\n")
})
