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
