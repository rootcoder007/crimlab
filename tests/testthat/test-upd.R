test_that("the version in the answer of the release page", {
  h <- c("HTTP/2 302 \r\n", "server: x\r\n",
         "Location: https://github.com/a/crimlab/releases/tag/v1.20.3\r\n")
  expect_identical(updver(h), "1.20.3")
  expect_identical(updver(sub("Location", "location", h, fixed = TRUE)),
                   "1.20.3")
  expect_identical(updver(c("HTTP/2 200 \r\n", "server: x\r\n")),
                   NA_character_)
  expect_identical(updver("location: https://x/releases/tag/v1.2;rm\r\n"),
                   NA_character_)
  expect_identical(updver(character()), NA_character_)
  # an address that cannot be reached is NA, without an error or a warning
  expect_no_warning(v <- updget("http://127.0.0.1:9/none"))
  expect_identical(v, NA_character_)
})

test_that("when the notice is off", {
  old <- Sys.getenv("CRIMLAB_NO_UPDATE_CHECK", NA)
  ask <- Sys.getenv("CRIMLAB_UPDATE_CHECK", NA)
  op <- options(crimlab.update = NULL)
  Sys.unsetenv("CRIMLAB_NO_UPDATE_CHECK")
  Sys.unsetenv("CRIMLAB_UPDATE_CHECK")
  # off unless asked for, interactive or not
  expect_true(updoff(FALSE))
  expect_true(updoff(TRUE))
  # asked for by the option, or by the variable
  options(crimlab.update = TRUE)
  expect_true(updoff(FALSE))
  expect_false(updoff(TRUE))
  options(crimlab.update = NULL)
  Sys.setenv(CRIMLAB_UPDATE_CHECK = "1")
  expect_false(updoff(TRUE))
  Sys.setenv(CRIMLAB_UPDATE_CHECK = "0")
  expect_true(updoff(TRUE))
  # FALSE and the older variable turn it off whatever else is set
  Sys.setenv(CRIMLAB_UPDATE_CHECK = "yes")
  options(crimlab.update = FALSE)
  expect_true(updoff(TRUE))
  options(crimlab.update = TRUE)
  Sys.setenv(CRIMLAB_NO_UPDATE_CHECK = "1")
  expect_true(updoff(TRUE))
  # tests are not interactive: nothing is asked, nothing is said
  expect_null(updrun())
  expect_no_message(.onAttach("", "crimlab"))
  options(op)
  if (is.na(old)) {
    Sys.unsetenv("CRIMLAB_NO_UPDATE_CHECK")
  } else {
    Sys.setenv(CRIMLAB_NO_UPDATE_CHECK = old)
  }
  if (is.na(ask)) {
    Sys.unsetenv("CRIMLAB_UPDATE_CHECK")
  } else {
    Sys.setenv(CRIMLAB_UPDATE_CHECK = ask)
  }
})

test_that("updchk() asks once and says what it found", {
  cur <- as.character(utils::packageVersion("crimlab"))
  local_mocked_bindings(updget = function(...) "999.0.0",
                        .package = "crimlab")
  expect_output(r <- updchk(), "crimlab 999.0.0 is available (this is",
                fixed = TRUE)
  expect_identical(r, list(current = cur, latest = "999.0.0", newer = TRUE))
  local_mocked_bindings(updget = function(...) cur, .package = "crimlab")
  expect_output(r <- updchk(), "is the latest release", fixed = TRUE)
  expect_false(r[["newer"]])
  local_mocked_bindings(updget = function(...) NA_character_,
                        .package = "crimlab")
  expect_output(r <- updchk(), "could not be reached", fixed = TRUE)
  expect_identical(r[["latest"]], NA_character_)
  expect_false(r[["newer"]])
  # an answer that is not a version is not shown
  local_mocked_bindings(updget = function(...) "1;x", .package = "crimlab")
  expect_output(r <- updchk(), "could not be reached", fixed = TRUE)
})

test_that("the notice, asked once a day and kept between", {
  d <- tempfile("upd")
  f <- file.path(d, "sub", "update")
  cnt <- new.env()
  cnt$n <- 0L
  get <- function() {
    cnt$n <- cnt$n + 1L
    "0.2.0"
  }
  m <- updmsg("0.1.9", f, now = 1000, get = get)
  expect_match(m, "^crimlab 0\\.2\\.0 is available \\(this is 0\\.1\\.9\\)")
  expect_match(m, "options(crimlab.update = FALSE) stops this notice",
               fixed = TRUE)
  expect_identical(readLines(f), "1000 0.2.0")
  # within the day the kept answer is used
  expect_identical(updmsg("0.1.9", f, now = 1000 + 86399, get = get), m)
  expect_identical(cnt$n, 1L)
  # the same or a later installed version: nothing
  expect_null(updmsg("0.2.0", f, now = 2000, get = get))
  expect_null(updmsg("0.10.0", f, now = 2000, get = get))
  # after a day, and with a clock set back, it is asked again
  expect_identical(updmsg("0.1.9", f, now = 1000 + 86400, get = get), m)
  expect_identical(updmsg("0.1.9", f, now = 5, get = get), m)
  expect_identical(cnt$n, 3L)
  # no answer is kept as none, and is not a notice
  expect_null(updmsg("0.1.9", f, now = 1e6, get = function() NA_character_))
  expect_identical(readLines(f), "1000000 NA")
  expect_null(updmsg("0.1.9", f, now = 1e6 + 1, get = get))
  # a kept file with anything else in it is asked again, never shown
  writeLines("x y", f)
  expect_identical(updmsg("0.1.9", f, now = 3e6, get = get), m)
  writeLines(paste(3e6, "9;system('x')"), f)
  expect_null(updmsg("0.1.9", f, now = 3e6 + 1, get = get))
  # a place that cannot be written: the answer is still used
  expect_identical(
    updmsg("0.1.9", file.path(f, "under-a-file", "update"), now = 1,
           get = get), m
  )
  unlink(d, recursive = TRUE)
})

test_that("the notice at attach", {
  local_mocked_bindings(updrun = function() "crimlab 9 is available",
                        .package = "crimlab")
  expect_message(.onAttach("", "crimlab"), "crimlab 9 is available",
                 fixed = TRUE, class = "packageStartupMessage")
})

test_that("the question goes to updget() unless a test gives its own", {
  local_mocked_bindings(updget = function(...) "0.3.0", .package = "crimlab")
  f <- tempfile("upd")
  expect_match(updmsg("0.1.0", f, now = 1e6),
               "crimlab 0.3.0 is available (this is 0.1.0)", fixed = TRUE)
  unlink(f)
})

test_that("the notice for this installation", {
  local_mocked_bindings(
    updoff = function(...) FALSE,
    updmsg = function(cur, file, ...) paste(cur, basename(file)),
    .package = "crimlab"
  )
  expect_identical(updrun(),
                   paste(utils::packageVersion("crimlab"), "update"))
})
