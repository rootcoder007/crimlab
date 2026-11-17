# The open data functions on a copy of the service in a temporary folder
# (helper-od.R): nothing here reaches the network.

test_that("odtab() loads a table as served and keeps it in the cache", {
  d <- odfix()
  odrun(d, {
    to <- getOption("timeout")
    x <- odtab("db1/t1")
    expect_identical(names(x), c("Case_ID", "offence type", "year"))
    expect_identical(x$year, c(2020L, 2021L))
    expect_identical(getOption("timeout"), to)
    kept <- file.path(tools::R_user_dir("crimlab", "cache"), "db1",
                      "t1.csv.gz")
    expect_true(file.exists(kept))
    # the kept copy matches the manifest: the service is not asked again
    served <- file.path(d, "db1", "t1.csv.gz")
    file.rename(served, paste0(served, ".bak"))
    expect_identical(odtab("db1/t1"), x)
    errs("odtab", list(
      list("could not be downloaded", "db1/t1", refresh = TRUE)
    ))
    expect_identical(getOption("timeout"), to)
    expect_true(file.exists(kept))
    # a kept copy that no longer matches is downloaded again
    file.rename(paste0(served, ".bak"), served)
    writeLines("x", kept)
    expect_identical(odtab("db1/t1"), x)
    expect_identical(.odhash(kept, "odtab"), .odhash(served, "odtab"))
  })
})

test_that("odtab() never uses a file that does not match the manifest", {
  d <- odfix()
  odgz(file.path(d, "db1", "t1.csv.gz"), c("Case_ID", "9"))
  odrun(d, {
    errs("odtab", list(list("does not match the SHA-256", "db1/t1")))
    kept <- file.path(tools::R_user_dir("crimlab", "cache"), "db1")
    expect_length(list.files(kept), 0L)
  })
})

test_that("odtab() matches columns without regard to case or underscores", {
  d <- odfix()
  odrun(d, {
    x <- odtab("db1/t1", cols = c("YEAR", "caseid", "offence type"))
    expect_identical(names(x), c("year", "Case_ID", "offence type"))
    expect_identical(x$Case_ID, 1:2)
    errs("odtab", list(
      list("no column `nope`", "db1/t1", cols = "nope"),
      list("`cols` must be column names", "db1/t1", cols = 1),
      list("`cols` must be column names", "db1/t1", cols = NA_character_)
    ))
  })
  # an exact name first; else one name only
  expect_identical(.odcol("ab", c("a_b", "ab"), "odtab"), 2L)
  expect_error(.odcol("AB", c("a_b", "Ab"), "odtab"),
               "^odtab: `AB` matches more than one column: a_b, Ab")
})

test_that("odtab() refuses in its own name", {
  d <- odfix()
  odrun(d, {
    errs("odtab", list(
      list("`db1/zz` is not a table", "db1/zz"),
      list("lists `db1/t2` without a SHA-256", "db1/t2"),
      list("`key` must be a single table key", 1),
      list("`key` must be a single table key", c("db1/t1", "db1/t2")),
      list("`refresh` must be TRUE or FALSE", "db1/t1", refresh = NA),
      list("`base` must be one", "db1/t1", base = "data.example"),
      list("manifest could not be downloaded", "db1/t1",
           base = odurl(file.path(d, "none")))
    ))
    odman(d, "<html></html>")
    errs("odtab", list(list("is not a manifest of tables", "db1/t1")))
    odman(d, r"({"datasets": system("x")})")
    errs("odtab", list(list("is not a manifest of tables", "db1/t1")))
  })
})

test_that("the manifest reader takes JSON and nothing else", {
  j <- .odjson(
    r"({"a": [1, -2.5e1, true, false, null], "b": "x\"\\\/yé", "c": {}})"
  )
  expect_identical(j$a, list(1, -25, TRUE, FALSE, NULL))
  expect_identical(j$b, "x\"\\/yé")
  expect_identical(j$c, list())
  bad <- c(r"({"a": system("x")})", "[e]", "[1] [2]", "[1,", "", r"(["\x"])")
  for (b in bad) expect_error(.odjson(b), info = b)
})

test_that("odls() lists the manifest's tables", {
  d <- odfix()
  gz <- file.path(d, "db1", "t1.csv.gz")
  odrun(d, {
    l <- odls()
    expect_identical(names(l), c("key", "rows", "bytes", "sha256", "columns"))
    expect_identical(l$key, c("db1/t1", "db1/t2"))
    expect_identical(l$rows, c(2, 1))
    expect_identical(l$bytes, c(file.size(gz), NA))
    expect_identical(l$sha256, c(.odhash(gz, "odls"), NA))
    expect_identical(
      l$columns,
      list(c("Case_ID", "offence type", "year"), character(0))
    )
    errs("odls", list(
      list("`base` must be one", base = NA),
      list("manifest could not be downloaded",
           base = odurl(file.path(d, "none")))
    ))
  })
})

test_that("odtab() refuses a table it cannot fetch, check, keep or read", {
  d <- odfix()
  gz <- file.path(d, "db1", "t1.csv.gz")
  mf <- file.path(d, "manifest.json")
  m <- readLines(mf)
  odrun(d, {
    odman(d, sub('"db": "db1", "table": "t1"', '"db": "..", "table": "t1"',
                 m, fixed = TRUE))
    errs("odtab", list(list("without a usable db and table", "db1/t1")))
    odman(d, m)
  })
  odrun(d, {
    odgz(gz, c("Case_ID,offence type,year", "1,theft,2020"))
    cur <- sub("[0-9a-f]{64}", .odhash(gz, "t"), m)
    odman(d, cur)
    expect_identical(nrow(odtab("db1/t1")), 1L)
    # a folder where the kept file goes: the download cannot replace it
    kept <- file.path(tools::R_user_dir("crimlab", "cache"), "db1",
                      "t1.csv.gz")
    unlink(kept)
    dir.create(file.path(kept, "x"), recursive = TRUE)
    suppressWarnings(errs("odtab", list(
      list("could not be kept", "db1/t1", refresh = TRUE)
    )))
    unlink(kept, recursive = TRUE)
    # an empty file: its SHA-256 is the manifest's, a table it is not
    odgz(gz, character(0))
    odman(d, sub("[0-9a-f]{64}", .odhash(gz, "t"), m))
    errs("odtab", list(list("`db1/t1` could not be read", "db1/t1")))
    # a cache folder that cannot be made
    blk <- file.path(d, "blk")
    writeLines("x", blk)
    Sys.setenv(R_USER_CACHE_DIR = blk)
    errs("odtab", list(list("cache folder .* cannot be made", "db1/t1")))
  })
})

test_that("a file too large to hash at once is refused", {
  skip_on_os("windows")
  skip_on_cran()
  f <- tempfile()
  on.exit(unlink(f))
  con <- file(f, "wb")
  seek(con, 2^31, rw = "write")
  writeBin(as.raw(1L), con)
  close(con)
  skip_if(is.na(file.size(f)) || file.size(f) <= .Machine$integer.max)
  expect_error(.odhash(f, "odtab"), "^odtab: .* is too large to check here")
})

test_that("a failed download is told in words", {
  u <- "https://data.example/db/t.csv.gz"
  w <- function(m) .odwhy(m, u, "`db/t`", 60)
  expect_match(w("cannot open URL: HTTP status was '404 Not Found'"),
               "does not serve `db/t` \\(404\\)")
  expect_match(w("HTTP status was '403 Forbidden'"), "refused `db/t` \\(403\\)")
  expect_match(w("HTTP status was '410 Gone'"), "\\(410\\): its publisher")
  expect_match(w("HTTP status was '429 Too Many Requests'"),
               "too many requests")
  expect_match(w("HTTP status was '503 Service Unavailable'"),
               "is unavailable \\(503\\)")
  expect_match(w("Forbidden"), "\\(403\\)")
  expect_match(w("Timeout was reached"), "did not answer within 60 seconds")
  expect_match(w("could not resolve host"),
               "could not be downloaded from https://data.example \\(could")
  expect_match(w(""), "\\(no reason given\\)")
  local_mocked_bindings(download.file = function(...) 1L, .package = "utils")
  errs("odls", list(list("no reason given", base = "https://data.example")))
})

test_that("a file that changed, with a manifest that lists no SHA-256 now", {
  d <- odfix()
  gz <- file.path(d, "db1", "t1.csv.gz")
  odrun(d, {
    st <- new.env()
    st$n <- 0L
    real <- .odsha
    local_mocked_bindings(.odsha = function(d) {
      st$n <- st$n + 1L
      if (st$n == 1L) real(d) else NA_character_
    })
    odgz(gz, c("a", "9"))
    errs("odtab", list(list("against none", "db1/t1")))
    expect_identical(st$n, 2L)
  })
})
