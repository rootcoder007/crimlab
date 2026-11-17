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
  expect_true(tryCatch(clprep(preset = bad), error = ok))
  expect_true(tryCatch(dastd(c(5, 3), c(100, 100), c("a", "b"),
                             ref = bad), error = ok))
  # the package's stop(), warning() and message(), and conditions given
  expect_true(tryCatch(stop("x ", bad), error = ok))
  expect_true(tryCatch(stop(simpleError(bad)), error = ok))
  expect_true(tryCatch(warning("x ", bad), warning = ok))
  expect_true(tryCatch(warning(simpleWarning(bad)), warning = ok))
  expect_true(tryCatch(message("x ", bad), message = ok))
  expect_true(tryCatch(message(simpleMessage(bad)), message = ok))
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

test_that(".clkey() keys rows apart unless they agree in every vector", {
  k <- .clkey(c("a\rb", "a", NA, "NA", "\\N", "a\\", "a"),
              c("c", "b\rc", "x", "x", "x", "\rc", "\\\rc"))
  expect_identical(anyDuplicated(k), 0L)
  expect_identical(.clkey(1:2, "x"), .clkey(c("1", "2"), "x"))
  expect_identical(.clkey(NA, NA), .clkey(NA_character_, NA_real_))
  expect_identical(.clkey(sep = "a", collapse = "b"), "a\rb")
})

test_that("the package's warning() and message() take base R's arguments", {
  f <- function() warning("imm", immediate. = TRUE, noBreaks. = TRUE)
  w <- tryCatch(f(), warning = identity)
  expect_identical(conditionMessage(w), "imm")
  expect_identical(conditionCall(w), quote(f()))
  # immediate. (warn 1 while it is signalled) leaves warn as it was
  op <- options(warn = 0)
  on.exit(options(op))
  suppressWarnings(f())
  expect_identical(as.integer(getOption("warn")), 0L)
  g <- function() message("m", appendLF = FALSE)
  m <- tryCatch(g(), message = identity)
  expect_identical(conditionMessage(m), "m")
  expect_identical(conditionCall(m), quote(message("m", appendLF = FALSE)))
  # base R's call is the message() call itself
  expect_identical(conditionCall(tryCatch(message("m"), message = identity)),
                   quote(message("m")))
})

test_that("a scalar keeps no dimensions; data frames and paths are checked", {
  expect_null(dim(.cknum(matrix(0.5), "p", "f")))
  expect_identical(.cktxtdf(data.frame(a = "x"), "d", "f"), data.frame(a = "x"))
  expect_error(.cktxtdf(data.frame(a = "\xff"), "d", "f"),
               "^f: `d` holds text that is not valid")
  expect_error(.cktxtdf(stats::setNames(data.frame(1), "\xfe"), "d", "f"),
               "^f: `d` holds text that is not valid")
  expect_error(.cktxtdf(data.frame(a = factor("\xfe")), "d", "f"),
               "^f: `d` holds text that is not valid")
  expect_null(.ckpath(NULL, "f"))
  expect_error(.ckpath("", "f"), "^f: `file` must be a file path, not empty")
  p <- file.path(tempdir(), "no", "such", "x.txt")
  expect_error(.ckpath(p, "f"), "^f: `file` is in a folder that does not")
  expect_identical(.ckpath(file.path(tempdir(), "x.txt"), "f"),
                   file.path(tempdir(), "x.txt"))
})
