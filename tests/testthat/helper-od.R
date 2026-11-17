# A copy of the open data service in a temporary folder, so that odtab()
# and odls() are tested without the network: a manifest and one table,
# read through a file:// address. odfix() makes the copy and gives its
# folder; odrun() runs code with the service at that folder and the
# cache in it, and restores both afterwards.
odfix <- function() {
  d <- tempfile("od")
  dir.create(file.path(d, "db1"), recursive = TRUE)
  gz <- file.path(d, "db1", "t1.csv.gz")
  odgz(gz, c("Case_ID,offence type,year", "1,theft,2020", "2,assault,2021"))
  odman(d, sprintf(
    r"({"generated_utc": "2026-01-01T00:00:00Z", "datasets": [
      {"db": "db1", "table": "t1", "key": "db1/t1", "rows": 2,
       "columns": ["Case_ID", "offence type", "year"], "bytes_gz": %.0f,
       "sha256": "%s", "source": "test", "meta": {"note": "café \/",
       "n": -1.5e2, "x": [true, false, null]}},
      {"db": "db1", "table": "t2", "key": "db1/t2", "rows": 1,
       "columns": [], "bytes_gz": null, "sha256": null, "source": "",
       "meta": {}}]})",
    file.size(gz),
    .odhash(gz, "odfix")
  ))
  d
}

odgz <- function(path, lines) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  writeLines(lines, con)
}

odman <- function(d, txt) writeLines(txt, file.path(d, "manifest.json"))

odurl <- function(d) {
  u <- normalizePath(d, winslash = "/", mustWork = FALSE)
  paste0("file://", if (!startsWith(u, "/")) "/", u)
}

odrun <- function(d, code) {
  old <- options(crimlab.data = odurl(d))
  cache <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
  Sys.setenv(R_USER_CACHE_DIR = file.path(d, "cache"))
  on.exit({
    options(old)
    if (is.na(cache)) {
      Sys.unsetenv("R_USER_CACHE_DIR")
    } else {
      Sys.setenv(R_USER_CACHE_DIR = cache)
    }
  })
  code
}
