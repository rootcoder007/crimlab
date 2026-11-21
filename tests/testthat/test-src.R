# The publisher's portals (R/src.R). Nothing here reaches the network: the
# download is replaced by one that serves the answers written below, in
# the shapes CKAN, Socrata, ArcGIS and Statistics Canada give.

# A download that serves `pages` (address -> text); any other address fails
# as a download does.
web <- function(pages, seen = NULL) {
  function(url, dest, wait, what, fn) {
    if (!is.null(seen)) {
      seen[["url"]] <- c(seen[["url"]], url)
    }
    if (is.null(pages[[url]])) {
      stop(fn, ": ", what, " could not be downloaded.", call. = FALSE)
    }
    writeLines(pages[[url]], dest)
    invisible(dest)
  }
}
serve <- function(pages, code, seen = NULL) {
  local_mocked_bindings(.odget = web(pages, seen), .package = "crimlab")
  code
}

scwds <- paste0(
  "https://www150.statcan.gc.ca/t1/wds/rest/",
  "getFullTableDownloadCSV/"
)
ck <- "https://portal.test"
ckpkg <- sprintf("%s/api/3/action/package_show?id=crimes", ck)
ckjson <- r"({"success": true, "result": {"name": "crimes", "resources": [
  {"name": "Readme", "format": "PDF", "url": "https://portal.test/r.pdf",
   "created": "2024-01-01", "last_modified": null, "size": null},
  {"name": "Crimes", "format": "csv", "url": "https://portal.test/c.csv",
   "created": "2024-01-01", "last_modified": "2026-02-03", "size": 1234},
  "not a resource"]}})"

test_that("ckanres() lists the files of a CKAN dataset", {
  r <- serve(
    stats::setNames(list(ckjson), ckpkg),
    ckanres(paste0(ck, "/"), "crimes")
  )
  expect_named(r, c("name", "format", "url", "modified", "bytes"))
  expect_identical(r[["name"]], c("Readme", "Crimes"))
  expect_identical(r[["url"]][2L], "https://portal.test/c.csv")
  expect_identical(r[["modified"]], c("2024-01-01", "2026-02-03"))
  expect_identical(r[["bytes"]], c(NA, 1234))
  # an answer that is not the dataset, or not JSON, or no answer
  no <- stats::setNames(
    list('{"success": false}', "<html>"),
    c(ckpkg, sprintf("%s/api/3/action/package_show?id=x%%20y", ck))
  )
  expect_error(
    serve(no, ckanres(ck, "crimes")),
    "ckanres: https://portal.test did not find the dataset `crimes`",
    fixed = TRUE
  )
  expect_error(
    serve(no, ckanres(ck, "x y")),
    paste(
      "ckanres: https://portal.test answered the dataset `x y`",
      "with something that is not JSON"
    ),
    fixed = TRUE
  )
  expect_error(serve(list(), ckanres(ck, "crimes")), "^ckanres: ")
  errs(
    "ckanres",
    list(
      list(
        "`portal` must be one address starting with http:// or https://",
        "portal.test",
        "a"
      ),
      list("`portal` must be one address", NA_character_, "a"),
      list("`id` must be one non-empty string", ck, "")
    )
  )
})

test_that("ckanls() searches a CKAN portal", {
  u <- sprintf("%s/api/3/action/package_search?q=use%%20of%%20force&rows=2", ck)
  j <- r"({"success": true, "result": {"count": 2, "results": [
    {"name": "uof", "title": "Use of force", "license_title": "OGL",
     "resources": [{}, {}]},
    {"name": "uof2", "title": "Use of force 2", "resources": []}]}})"
  r <- serve(stats::setNames(list(j), u), ckanls(ck, "use of force", rows = 2))
  expect_identical(r[["id"]], c("uof", "uof2"))
  expect_identical(r[["title"]], c("Use of force", "Use of force 2"))
  expect_identical(r[["files"]], c(2L, 0L))
  expect_identical(r[["licence"]], c("OGL", NA))
  errs(
    "ckanls",
    list(
      list("`q` must be one non-empty string", ck, 1),
      list("`rows` must be a whole number from 1 to 1000", ck, "a", rows = 0),
      list("`rows` must be a whole number from 1 to 1000", ck, "a", rows = 2.5),
      list("`rows` must be a whole number from 1 to 1000", ck, "a", rows = 1001)
    )
  )
})

sj <- r"({"id": "ijzp-q8t2", "name": "Crimes", "rowsUpdatedAt": 1700000000,
  "columns": [{"name": "ID"}, {"name": "Date"}, "x"]})"

test_that("socres() describes a Socrata dataset", {
  u <- "https://data.city.test/api/views/ijzp-q8t2.json"
  r <- serve(
    stats::setNames(list(sj), u),
    socres("data.city.test", "ijzp-q8t2")
  )
  expect_identical(r[["name"]], "Crimes")
  expect_identical(
    r[["url"]],
    "https://data.city.test/api/views/ijzp-q8t2/rows.csv?accessType=DOWNLOAD"
  )
  expect_identical(r[["columns"]], c("ID", "Date"))
  expect_identical(as.numeric(r[["updated"]]), 1700000000)
  expect_identical(attr(r[["updated"]], "tzone"), "UTC")
  # with the scheme given, and an answer without a date
  nod <- '{"id": "ijzp-q8t2", "name": "Crimes", "columns": []}'
  r2 <- serve(
    stats::setNames(list(nod), u),
    socres("https://data.city.test", "ijzp-q8t2")
  )
  expect_true(is.na(r2[["updated"]]))
  expect_identical(r2[["columns"]], character(0))
  expect_error(
    serve(
      stats::setNames(list('{"id": "aaaa-bbbb"}'), u),
      socres("data.city.test", "ijzp-q8t2")
    ),
    "socres: https://data.city.test did not describe the dataset `ijzp-q8t2`",
    fixed = TRUE
  )
  expect_error(
    serve(
      stats::setNames(list("[1]"), u),
      socres("data.city.test", "ijzp-q8t2")
    ),
    "did not describe the dataset",
    fixed = TRUE
  )
  errs(
    "socres",
    list(
      list("`domain` must be one address", 1, "ijzp-q8t2"),
      list("`id` must be a Socrata identifier", "data.city.test", "crimes"),
      list("`id` must be one non-empty string", "data.city.test", NA_character_)
    )
  )
})

lay <- "https://gis.test/rest/services/Crime/FeatureServer/0"
aj <- r"({"name": "Major crime", "maxRecordCount": 2000,
  "editingInfo": {"lastEditDate": 1700000000000},
  "fields": [{"name": "OBJECTID"}, {"name": "MCI"}]})"

test_that("arcres() describes an ArcGIS layer", {
  p <- stats::setNames(
    list(aj, '{"count": 420000}'),
    c(
      paste0(lay, "?f=json"),
      sprintf("%s/query?where=1%%3D1&returnCountOnly=true&f=json", lay)
    )
  )
  r <- serve(p, arcres(lay))
  expect_identical(r[["name"]], "Major crime")
  expect_identical(
    r[["url"]],
    sprintf("%s/query?where=1%%3D1&outFields=*&f=geojson", lay)
  )
  expect_identical(r[["columns"]], c("OBJECTID", "MCI"))
  expect_identical(as.numeric(r[["updated"]]), 1700000000)
  expect_identical(r[["rows"]], 420000)
  expect_identical(r[["max"]], 2000)
  # no edit date, a count that is not one, no limit
  p2 <- stats::setNames(
    list('{"name": "L", "fields": []}', "[]"),
    names(p)
  )
  r2 <- serve(p2, arcres(lay))
  expect_true(is.na(r2[["updated"]]))
  expect_true(is.na(r2[["rows"]]))
  expect_true(is.na(r2[["max"]]))
  bad <- stats::setNames(list('{"error": {"code": 400}}'), names(p)[1L])
  expect_error(
    serve(bad, arcres(lay)),
    paste0("arcres: ", lay, " did not describe a layer"),
    fixed = TRUE
  )
  errs(
    "arcres",
    list(
      list("`layer` must be one address", "gis.test/FeatureServer/0"),
      list(
        "`layer` must end in a layer number",
        "https://gis.test/Crime/FeatureServer"
      )
    )
  )
})

test_that("scres() asks Statistics Canada where a table is", {
  u <- paste0(scwds, "35100177/en")
  ok <- paste0(
    '{"status": "SUCCESS", "object": ',
    '"https://www150.statcan.gc.ca/n1/tbl/csv/35100177-eng.zip"}'
  )
  p <- stats::setNames(list(ok), u)
  want <- "https://www150.statcan.gc.ca/n1/tbl/csv/35100177-eng.zip"
  expect_identical(serve(p, scres("35100177")), want)
  expect_identical(serve(p, scres("35-10-0177-01")), want)
  expect_identical(serve(p, scres(35100177)), want)
  fr <- sub("/en$", "/fr", u)
  expect_identical(
    serve(stats::setNames(list(ok), fr), scres("3510017701", lang = "fr")),
    want
  )
  for (j in c(
    '{"status": "FAILED", "object": "x"}',
    "[1]",
    '{"status": "SUCCESS", "object": "ftp://x"}',
    '{"status": "SUCCESS"}'
  )) {
    expect_error(
      serve(stats::setNames(list(j), u), scres("35100177")),
      "scres: Statistics Canada did not find the table `35100177`",
      fixed = TRUE
    )
  }
  errs(
    "scres",
    list(
      list("`pid` must be a table number of eight digits", "3510"),
      list("`pid` must be a table number of eight digits", "abcdefgh"),
      list("`pid` must be one non-empty string", NA),
      list("`lang` must be one of", "35100177", lang = "de")
    )
  )
})

test_that("srcls() is the list of sources, each key once", {
  s <- srcls()
  expect_named(s, c("key", "publisher", "what", "kind", "portal", "id"))
  expect_false(anyDuplicated(s[["key"]]) > 0L)
  expect_true(all(vapply(s, is.character, NA)))
  expect_true(all(
    nzchar(s[["key"]]) &
      nzchar(s[["publisher"]]) &
      nzchar(s[["what"]]) &
      nzchar(s[["portal"]])
  ))
  expect_setequal(
    unique(s[["kind"]]),
    c("arcgis", "ckan", "page", "socrata", "statcan", "url")
  )
  # a page and a file are an address; every other kind names a dataset or
  # is a portal to search
  expect_false(any(nzchar(s[["id"]][s[["kind"]] %in% c("page", "url")])))
  expect_true(all(grepl(
    "^[a-z0-9]{4}-[a-z0-9]{4}$",
    s[["id"]][s[["kind"]] == "socrata"]
  )))
  expect_true(all(grepl("^[0-9]{8}$", s[["id"]][s[["kind"]] == "statcan"])))
  expect_true(all(
    c(
      "tps-mci",
      "statcan-crime",
      "nypd-complaints",
      "cpd-crimes",
      "uk-police",
      "eurostat-crime",
      "portal-australia"
    ) %in%
      s[["key"]]
  ))
})

test_that("srcurl() resolves a key at the publisher's portal", {
  s <- srcls()
  at <- function(k) s[s[["key"]] == k, ]
  # a page, a file and a portal: the address itself
  expect_identical(srcurl("siu-reports"), at("siu-reports")[["portal"]])
  expect_identical(
    srcurl("uk-police"),
    "https://data.police.uk/data/archive/latest.zip"
  )
  expect_identical(srcurl("portal-ontario"), "https://data.ontario.ca")
  # Statistics Canada
  su <- paste0(scwds, "35100177/en")
  expect_identical(
    serve(
      stats::setNames(
        list('{"status": "SUCCESS", "object": "https://x.test/t.zip"}'),
        su
      ),
      srcurl("statcan-crime")
    ),
    "https://x.test/t.zip"
  )
  # CKAN: the CSV when there is one, else the first file
  o <- at("ontario-uof")
  cu <- paste0(o[["portal"]], "/api/3/action/package_show?id=", o[["id"]])
  expect_identical(
    serve(stats::setNames(list(ckjson), cu), srcurl("ontario-uof")),
    "https://portal.test/c.csv"
  )
  pdf <- sub('"format": "csv"', '"format": "XLSX"', ckjson, fixed = TRUE)
  expect_identical(
    serve(stats::setNames(list(pdf), cu), srcurl("ontario-uof")),
    "https://portal.test/r.pdf"
  )
  none <- '{"success": true, "result": {"resources": []}}'
  expect_error(
    serve(stats::setNames(list(none), cu), srcurl("ontario-uof")),
    paste(
      "srcurl: https://data.ontario.ca lists no file for",
      "`use-of-force-in-correctional-institutions`"
    ),
    fixed = TRUE
  )
  # a refusal of the portal comes under srcurl()'s name
  expect_error(serve(list(), srcurl("ontario-uof")), "^srcurl: ")
  # Socrata
  cp <- at("cpd-crimes")
  expect_identical(
    serve(
      stats::setNames(
        list(sj),
        "https://data.cityofchicago.org/api/views/ijzp-q8t2.json"
      ),
      srcurl("cpd-crimes")
    ),
    paste0(
      "https://data.cityofchicago.org/api/views/ijzp-q8t2/",
      "rows.csv?accessType=DOWNLOAD"
    )
  )
  expect_identical(cp[["id"]], "ijzp-q8t2")
  # ArcGIS
  t <- at("tps-mci")
  l <- sprintf("%s/%s/FeatureServer/0", t[["portal"]], t[["id"]])
  p <- stats::setNames(
    list(aj, '{"count": 1}'),
    c(
      paste0(l, "?f=json"),
      sprintf("%s/query?where=1%%3D1&returnCountOnly=true&f=json", l)
    )
  )
  expect_identical(
    serve(p, srcurl("tps-mci")),
    sprintf("%s/query?where=1%%3D1&outFields=*&f=geojson", l)
  )
  errs(
    "srcurl",
    list(
      list("`nokey` is not a key of srcls()", "nokey"),
      list("`key` must be one non-empty string", 1)
    )
  )
})

test_that("srcget() downloads to a file and says what it wrote", {
  d <- tempfile("srcget")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  src <- file.path(d, "a.csv")
  writeLines(c("id,n", "1,2"), src)
  u <- paste0(
    "file://",
    if (!startsWith(src, "/")) "/",
    normalizePath(src, winslash = "/")
  )
  out <- file.path(d, "b.csv")
  r <- srcget(u, out)
  expect_named(r, c("path", "bytes", "sha256", "url"))
  expect_identical(readLines(out), c("id,n", "1,2"))
  expect_identical(r[["path"]], out)
  expect_identical(r[["bytes"]], file.size(src))
  expect_identical(r[["sha256"]], .odhash(src, "srcget"))
  expect_identical(r[["url"]], u)
  expect_invisible(srcget(u, out, over = TRUE))
  expect_false(file.exists(paste0(out, ".part")))
  # an existing file stays unless told otherwise; a failure leaves no part
  expect_error(
    srcget(u, out),
    "exists; `over = TRUE` replaces it",
    fixed = TRUE
  )
  expect_error(srcget(paste0(u, ".none"), file.path(d, "c.csv")), "^srcget: ")
  expect_false(file.exists(file.path(d, "c.csv.part")))
  # a download that leaves nothing to move into place
  local_mocked_bindings(
    .odget = function(url, dest, ...) {
      writeLines("x", file.path(d, "elsewhere"))
    },
    .package = "crimlab"
  )
  expect_error(
    srcget("https://x.test/f", file.path(d, "g")),
    "could not be written",
    fixed = TRUE
  )
  errs(
    "srcget",
    list(
      list("`url` must be one address", "x.test/f", out),
      list("`over` must be TRUE or FALSE", u, out, over = NA),
      list("`dest` must be one file path", u, ""),
      list("`dest` must be one file path", u, c("a", "b")),
      list("`wait` must be one positive number of seconds", u, out, wait = 0),
      list("`wait` must be one positive number of seconds", u, out, wait = "a")
    )
  )
})

test_that("ovrget() saves a director's report by its number", {
  d <- tempfile("ovrget")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  seen <- new.env()
  en <- "https://www.siu.on.ca/en/directors_report_details.php?drid=4100"
  fr <- "https://www.siu.on.ca/fr/directors_report_details.php?drid=4100"
  p <- stats::setNames(
    list("<html>report</html>", "<html>rapport</html>"),
    c(en, fr)
  )
  r <- serve(p, ovrget(4100, file.path(d, "r.html")), seen)
  expect_identical(r[["url"]], en)
  expect_identical(readLines(file.path(d, "r.html")), "<html>report</html>")
  serve(p, ovrget(4100, file.path(d, "f.html"), lang = "fr"), seen)
  expect_identical(seen[["url"]], c(en, fr))
  expect_error(
    serve(p, ovrget(4100, file.path(d, "r.html"))),
    "ovrget: .*exists; `over = TRUE` replaces it"
  )
  errs(
    "ovrget",
    list(
      list("`drid` must be one positive report number", 0, "x"),
      list("`drid` must be one positive report number", 1.5, "x"),
      list("`drid` must be one positive report number", "4100", "x"),
      list("`lang` must be one of", 1, "x", lang = "de")
    )
  )
})

# The Indigenous tables of the open data service, on a copy of it.
indfix <- function() {
  d <- tempfile("ind")
  dir.create(file.path(d, "statcan_indigenous"), recursive = TRUE)
  gz <- file.path(d, "statcan_indigenous", "t1.csv.gz")
  odgz(gz, c("geo,value", "Canada,1"))
  odman(
    d,
    sprintf(
      r"({"generated_utc": "2026-01-01T00:00:00Z", "datasets": [
      {"db": "statcan_indigenous", "table": "t1",
       "key": "statcan_indigenous/t1",
       "rows": 1, "columns": ["geo", "value"], "bytes_gz": %.0f, "sha256": "%s",
       "source": "statcan", "meta": {"description": "Perceived mental health",
       "license": "OGL-Canada", "source": "Statistics Canada"}},
      {"db": "statcan_indigenous", "table": "t2",
       "key": "statcan_indigenous/t2",
       "rows": 5, "columns": [], "bytes_gz": null, "sha256": null,
       "source": "statcan", "meta": {"description": "Population by identity"}},
      {"db": "census2021_ipp", "table": "t3", "key": "census2021_ipp/t3",
       "rows": 9, "columns": [], "bytes_gz": null, "sha256": null,
       "source": "statcan", "meta": "none"},
      {"db": "census2021_ipp", "table": "t4", "key": "census2021_ipp/t4",
       "rows": 2, "columns": [], "bytes_gz": null, "sha256": null,
       "source": "statcan",
       "meta": {"description": "Population and dwellings"}},
      {"db": "otis", "table": "a01", "key": "otis/a01", "rows": 3,
       "columns": [], "bytes_gz": null, "sha256": null, "source": "x",
       "meta": {"description": "Population in custody"}},
      "not a dataset"]})",
      file.size(gz),
      .odhash(gz, "indfix")
    )
  )
  d
}

test_that("indtabs() lists the Indigenous tables of the service", {
  d <- indfix()
  odrun(d, {
    t <- indtabs()
    expect_named(t, c("key", "rows", "description", "licence", "source"))
    expect_identical(
      t[["key"]],
      c(
        "statcan_indigenous/t1",
        "statcan_indigenous/t2",
        "census2021_ipp/t3",
        "census2021_ipp/t4"
      )
    )
    expect_identical(t[["rows"]], c(1, 5, 9, 2))
    expect_identical(
      t[["description"]],
      c(
        "Perceived mental health",
        "Population by identity",
        "",
        "Population and dwellings"
      )
    )
    expect_identical(t[["licence"]], c("OGL-Canada", "", "", ""))
    expect_identical(indtabs("MENTAL")[["key"]], "statcan_indigenous/t1")
    expect_identical(nrow(indtabs("no such thing")), 0L)
    expect_identical(indtabs(base = odurl(d)), t)
  })
  errs(
    "indtabs",
    list(
      list("`pat` must be one non-empty string", ""),
      list("`pat` is not a valid regular expression", "("),
      list(
        "`base` must be one http://, https:// or file:// address",
        NULL,
        base = "x"
      )
    )
  )
})

test_that("indiex() finds the tables on a subject and loads one", {
  d <- indfix()
  odrun(d, {
    expect_identical(indiex("mentalhealth")[["key"]], "statcan_indigenous/t1")
    # the census subjects are the census tables alone
    expect_identical(indiex("censuspop")[["key"]], "census2021_ipp/t4")
    expect_identical(nrow(indiex("suicide")), 0L)
    x <- indiex("mentalhealth", load = TRUE)
    expect_identical(x, odtab("statcan_indigenous/t1"))
    expect_named(x, c("geo", "value"))
    expect_error(
      indiex("suicide", load = TRUE),
      "indiex: the service holds no table on `suicide`",
      fixed = TRUE
    )
    # a refusal of the service comes under indiex()'s name
    expect_error(indiex("censuspop", load = TRUE), "^indiex: ")
    expect_identical(indiex()[["key"]], character(0))
  })
  errs(
    "indiex",
    list(
      list("`name` must be one of", "nothing"),
      list("`load` must be TRUE or FALSE", "health", load = NA)
    )
  )
})
