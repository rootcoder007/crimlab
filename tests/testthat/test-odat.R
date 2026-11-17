# The open data functions on a copy of the service in a temporary folder
# (helper-od.R): nothing here reaches the network.

test_that("odtab() loads a table as served and keeps it in the cache", {
  d <- odfix()
  odrun(d, {
    to <- getOption("timeout")
    x <- odtab("db1/t1")
    expect_named(x, c("Case_ID", "offence type", "year"))
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
    expect_named(x, c("year", "Case_ID", "offence type"))
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
  expect_identical(j$c, stats::setNames(list(), character()))
  bad <- c(r"({"a": system("x")})", "[e]", "[1] [2]", "[1,", "", r"(["\x"])")
  for (b in bad) expect_error(.odjson(b), info = b)
})

test_that("odls() lists the manifest's tables", {
  d <- odfix()
  gz <- file.path(d, "db1", "t1.csv.gz")
  odrun(d, {
    l <- odls()
    expect_named(l, c("key", "rows", "bytes", "sha256", "columns"))
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
    cur <- sub(r"("rows": 2)", r"("rows": 1)", cur, fixed = TRUE)
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
    errs("odtab", list(list("`db1/t1` is empty: it has no header line",
                            "db1/t1")))
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

# A mirror holding one table db1/t1 of the given CSV lines, its manifest
# listing rows (the lines after the header unless given).
odone <- function(lines, rows = length(lines) - 1L) {
  d <- tempfile("od")
  dir.create(file.path(d, "db1"), recursive = TRUE)
  gz <- file.path(d, "db1", "t1.csv.gz")
  odgz(gz, lines)
  odman(d, sprintf(
    r"({"datasets": [{"db": "db1", "table": "t1", "key": "db1/t1",
      "rows": %s, "sha256": "%s"}]})",
    rows, .odhash(gz, "odone")
  ))
  d
}

test_that("odtab() gives the values the file holds", {
  d <- odone(c(
    "id,sex,x,n,day,big,e,neg,miss,txt",
    "007,F,1.5,2020,2024-01-31,12345678901234567890,1e400,-0.25,,NA",
    "7,F,2,2021,2024-02-01,1,2,-1,,a",
    "0007,F,10,-3,2024-02-02,2,3,0,,b"
  ))
  odrun(d, {
    x <- odtab("db1/t1")
    expect_identical(x$id, c("007", "7", "0007"))
    expect_identical(x$sex, c("F", "F", "F"))
    expect_identical(x$x, c(1.5, 2, 10))
    expect_identical(x$n, c(2020L, 2021L, -3L))
    expect_identical(x$day, c("2024-01-31", "2024-02-01", "2024-02-02"))
    expect_identical(x$big, c("12345678901234567890", "1", "2"))
    expect_identical(x$e, c("1e400", "2", "3"))
    expect_identical(x$neg, c(-0.25, -1, 0))
    expect_identical(x$miss, rep(NA_character_, 3))
    expect_identical(x$txt, c(NA, "a", "b"))
    expect_identical(odtab("db1/t1", cols = c("n", "id")),
                     x[c("n", "id")])
  })
  d <- odone(c("w,v", "1,3000000000", "2,1.0", "3,2.5e2"))
  odrun(d, {
    x <- odtab("db1/t1")
    expect_identical(x$w, 1:3)
    expect_identical(x$v, c(3e9, 1, 250))
  })
})

test_that("odtab() refuses a ragged file and a count not the manifest's", {
  long <- c("a,b", "1,2", "3,4", "5,6", "7,8", "9,10", "11,12,13,14")
  odrun(odone(long), errs("odtab", list(
    list("line 7 of `db1/t1` has 4 fields where its header has 2", "db1/t1")
  )))
  odrun(odone(c("a,b", "1,2", "3")), errs("odtab", list(
    list("line 3 of `db1/t1` has 1 fields where its header has 2", "db1/t1")
  )))
  odrun(odone(c("a,b", "1,2", "3,4"), rows = 3), errs("odtab", list(
    list("`db1/t1` has 2 rows where its manifest lists 3", "db1/t1")
  )))
  # a field in quotes over two lines, and a blank line, are one row each
  d <- odone(c("a,b", "\"x", "y\",2", "", "z,3"), rows = 2)
  odrun(d, expect_identical(odtab("db1/t1")$a, c("x\ny", "z")))
})

test_that("the manifest reader takes all JSON, of any length", {
  j <- .odjson(r"({"": 1, "a": {"": [2]}})")
  expect_named(j, c("", "a"))
  expect_identical(j[[1L]], 1)
  expect_identical(j$a, stats::setNames(list(list(2)), ""))
  expect_identical(.odjson(r"([{"": [null, {}]}])"),
                   list(stats::setNames(list(list(
                     NULL, stats::setNames(list(), character())
                   )), "")))
  s <- strrep("ab", 5e6)
  expect_identical(.odjson(paste0(r"({"s": ")", s, '"}'))$s, s)
  e <- strrep(r"(\n)", 1e5)
  expect_identical(nchar(.odjson(paste0('["', e, '"]'))[[1L]]), 100000L)
  expect_identical(.odjson(r"(["\ud83d\ude00", "\\ud83d"])"),
                   list("\U1F600", "\\ud83d"))
  # a lone half of a surrogate pair, U+0000 (no R text holds it) and bytes
  # that are not UTF-8 are not text: refused
  bad <- c(r"({"a"})", r"({"a", 1})", r"(["a": 1])", r"({"a": "b": 1})",
           r"("a": 1)", r"(["\q"])", "[tru]", "[1e]",
           r"({"": 1, 2})", r"(["": 1])", r"({"": "": 1})", r"("": 1)",
           r"(["\ud83d x"])", r"(["\u0000"])", "[\"\xff\"]")
  for (b in bad) expect_error(.odjson(b), info = b)
})

# a table of the bytes b, gzip'd, listed with rows rows
odraw <- function(b, rows) {
  d <- tempfile("od")
  dir.create(file.path(d, "db1"), recursive = TRUE)
  gz <- file.path(d, "db1", "t1.csv.gz")
  con <- gzfile(gz, "wb")
  writeBin(b, con)
  close(con)
  odman(d, sprintf(
    r"({"datasets": [{"db": "db1", "table": "t1", "key": "db1/t1",
      "rows": %s, "sha256": "%s"}]})",
    rows, .odhash(gz, "odraw")
  ))
  d
}

test_that("odtab() reads the CSV it is given, row for row, names as served", {
  # a blank line is a missing value of a one-column table
  odrun(odone(c("v", "1", "", "3")),
        expect_identical(odtab("db1/t1")$v, c(1L, NA, 3L)))
  odrun(odone(c("a,a", "1,2")),
        expect_named(odtab("db1/t1"), c("a", "a")))
  odrun(odone(c(",x,", "1,2,3")),
        expect_named(odtab("db1/t1"), c("", "x", "")))
  # numbers only when held as written: underflow stays text
  d <- odone(c("u,s,n", "1e-400,1e-320,1e-300", "1,1,2"))
  odrun(d, {
    x <- odtab("db1/t1")
    expect_identical(x$u, c("1e-400", "1"))
    expect_identical(x$s, c("1e-320", "1"))
    expect_identical(x$n, c(1e-300, 2))
  })
  # quotes, a byte order mark, CRLF; bytes that are not UTF-8 as <xx>
  b <- c(as.raw(c(0xef, 0xbb, 0xbf)),
         charToRaw("a,b\r\n\"x,\"\"y\"\"\nz\",2\r\n"),
         charToRaw("caf"), as.raw(0xe9), charToRaw(",3\r\n"))
  odrun(odraw(b, 2), {
    x <- expect_silent(odtab("db1/t1"))
    expect_identical(x$a, c("x,\"y\"\nz", "caf<e9>"))
    expect_identical(x$b, 2:3)
  })
})

test_that("odtab() refuses what is not CSV, naming the line", {
  bad <- list(
    list(c("a,b", "1,x\"y"), "line 2 of `db1/t1` has a quote inside"),
    list(c("a,b", "1,2", "\"3,4"), "line 3 of `db1/t1` opens a quoted"),
    list(c("a", "\"x\"y"), "line 2 of `db1/t1` has text after the quote")
  )
  for (k in bad) {
    expect_no_warning(odrun(odone(k[[1L]], rows = 1),
                            errs("odtab", list(list(k[[2L]], "db1/t1")))))
  }
  odrun(odraw(c(charToRaw("a\n1"), as.raw(0), charToRaw("\n")), 1),
        errs("odtab", list(list("line 2 of `db1/t1` holds a NUL", "db1/t1"))))
  odrun(odraw(raw(), 0),
        errs("odtab", list(list("`db1/t1` is empty", "db1/t1"))))
})

test_that("the manifest reader is JSON's grammar, of any depth", {
  # every case of RFC 8259's grammar the round-2 conformance set names
  # valid parses, every invalid one is refused
  ok <- c("{}", "[]", "[-0]", "[1E5]", "[-1.25e-3]", r"(["\/\b\f\n\r\t"])",
          r"({"":{"":[]}})", " \t\r\n[ 1 , 2 ]\n ", r"("str")", "1", "null",
          r"({"a":1,"a":2})", "[1.7976931348623157e308]", "[5e-324]")
  for (s in ok) expect_no_error(.odjson(s))
  no <- c("[01]", "[1.]", "[.5]", "[+1]", "[1,]", "[1 2]", "{a:1}", "[NaN]",
          "[1e]", "[3-1]", "[3 -1]", "[}", "{]", "[", "", "[\"a\tb\"]",
          r"(["\x41"])", r"(["\u12"])", "[1]x", "{}{}", "[1]\n[2]", "﻿[1]",
          "[1]\001", r"(["\ud800"])")
  for (s in no) expect_error(.odjson(s), "not JSON", info = s)
  expect_identical(.odjson("[1e400]"), list(Inf))
  # 10000 levels inside the outermost are read, 10001 refused; nothing
  # is evaluated
  d <- paste0(strrep("[", 10002L), strrep("]", 10002L))
  expect_length(.odjson(substr(d, 2L, nchar(d) - 1L)), 1L)
  expect_error(.odjson(d), "not JSON")
  expect_identical(.odjson(r"-(["system('x')", "`a`"])-"),
                   list("system('x')", "`a`"))
})
