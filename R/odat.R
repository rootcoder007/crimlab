# The tables of the project's open data service. The service's manifest,
# read at every call, says all there is to know about a table: its key,
# rows, columns, size and the SHA-256 of its gzip file as served. Nothing
# about any table is written here, so a publisher's new file needs no new
# release. A file is used only when its SHA-256 is the manifest's.

#' Tables of the open data service
#'
#' `odtab()` loads a table of the project's open data service by its key,
#' as the service serves it at the time of the call. Everything about a
#' table (its address, rows, columns and SHA-256) is read from the
#' service's manifest at run time, so a table that its publisher has
#' replaced is loaded as it is now, with no new release of crimlab.
#'
#' A downloaded file is used only when its SHA-256 is the one the
#' manifest lists. When it is not, the manifest is read once more (the
#' table may have been replaced between the two requests) and, if the
#' file still does not match, it is deleted and `odtab()` stops. A table
#' the manifest lists without a SHA-256 cannot be checked and is not
#' loaded.
#'
#' Tables are kept in `tools::R_user_dir("crimlab", "cache")`. A kept
#' file is used again only while its SHA-256 is the one the manifest lists
#' now; otherwise, or with `refresh = TRUE`, the table is downloaded again.
#'
#' The columns are those of the file, named as served. `cols` picks some
#' of them, matched without regard to case or underscores (`"case_id"`
#' finds `"CaseID"`), so that a publisher's respelling of a name does not
#' break a call; an exact name is preferred.
#'
#' Errors: a key the manifest does not list; a table listed without a
#' SHA-256; a file that does not match it; a column the table does not
#' have, or a name that matches more than one; a service that cannot be
#' reached, answers with an error or does not answer in time (60 seconds
#' for the manifest, 600 for a table, or `getOption("timeout")` if
#' longer).
#'
#' @param key A table's key, as the service's manifest lists it.
#' @param cols The names of the columns to keep, or `NULL` (all).
#' @param refresh `TRUE` downloads the table again even when the cache
#'   holds a copy that matches.
#' @param base The address of the service, or of a mirror of it (a folder
#'   holding `manifest.json` and `<db>/<table>.csv.gz`, which may be a
#'   `file://` address). The default is `getOption("crimlab.data")`, or
#'   the project's service when that is not set.
#' @return The table as a data frame.
#' @export
odtab <- function(
  key,
  cols = NULL,
  refresh = FALSE,
  base = getOption("crimlab.data", "https://data.rmorie.com")
) {
  fn <- "odtab"
  if (!is.character(key) || length(key) != 1L || is.na(key)) {
    .ckstop(fn, "`key` must be a single table key")
  }
  refresh <- .ckflag(refresh, "refresh", fn)
  base <- .odbase(base, fn)
  d <- .odent(.odman(base, fn), key, base, fn)
  want <- .odsha(d)
  if (is.na(want)) {
    .ckstop(
      fn,
      "%s lists `%s` without a SHA-256: it cannot be checked, so is not used",
      base,
      key
    )
  }
  dir <- file.path(tools::R_user_dir("crimlab", "cache"), d[["db"]])
  path <- file.path(dir, paste0(d[["table"]], ".csv.gz"))
  if (refresh || !file.exists(path) || .odhash(path, fn) != want) {
    .odfetch(d, want, base, dir, path, fn)
  }
  .odread(path, cols, key, fn)
}

# The service's address: one http(s):// or file:// address, without the
# trailing slashes.
.odbase <- function(base, fn) {
  ok <- is.character(base) && length(base) == 1L && !is.na(base) &&
    grepl("^(https?|file)://.", base)
  if (!ok) {
    .ckstop(
      fn,
      "`base` must be one http://, https:// or file:// address"
    )
  }
  sub("/+$", "", base)
}

# The manifest, as lists.
.odman <- function(base, fn) {
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f), add = TRUE)
  .odget(paste0(base, "/manifest.json"), f, 60, "the manifest", fn)
  m <- tryCatch(
    .odjson(paste(readLines(f, warn = FALSE, encoding = "UTF-8"),
                  collapse = "\n")),
    error = function(e) NULL
  )
  if (!is.list(m) || !is.list(m[["datasets"]])) {
    .ckstop(fn, "%s/manifest.json is not a manifest of tables", base)
  }
  m
}

# The manifest's entry for key, with a db and table that are safe as a
# file name and in an address.
.odent <- function(m, key, base, fn) {
  for (d in m[["datasets"]]) {
    if (is.list(d) && identical(d[["key"]], key)) {
      ok <- vapply(d[c("db", "table")], function(x) {
        is.character(x) && length(x) == 1L &&
          grepl("^[A-Za-z0-9_.-]+$", x) && !x %in% c(".", "..")
      }, NA)
      if (!all(ok)) {
        .ckstop(fn, "%s lists `%s` without a usable db and table", base, key)
      }
      return(d)
    }
  }
  .ckstop(fn, "`%s` is not a table of %s", key, base)
}

# The entry's SHA-256 in lower case, or NA when it lists none.
.odsha <- function(d) {
  s <- d[["sha256"]]
  if (!is.character(s) || length(s) != 1L || is.na(s)) {
    return(NA_character_)
  }
  s <- sub("^sha256:", "", tolower(s))
  if (grepl("^[0-9a-f]{64}$", s)) s else NA_character_
}

# The SHA-256 of a file, in lower-case hex.
# ponytail: the whole file is read at once, so up to 2^31 - 1 bytes; hash
# in chunks if a table's gzip file ever outgrows that.
.odhash <- function(path, fn) {
  n <- file.size(path)
  if (is.na(n) || n > .Machine$integer.max) {
    .ckstop(fn, "%s is too large to check here", basename(path))
  }
  paste(.Call(C_clab_r_sha256, readBin(path, "raw", n)), collapse = "")
}

# The table downloaded into the cache. A file whose SHA-256 is not the
# manifest's is checked against the manifest read once more, then
# deleted and refused; only a match replaces the kept copy.
.odfetch <- function(d, want, base, dir, path, fn) {
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE,
                                      showWarnings = FALSE)) {
    .ckstop(fn, "the cache folder %s cannot be made", dir)
  }
  tmp <- tempfile(tmpdir = dir, fileext = ".part")
  on.exit(unlink(tmp), add = TRUE)
  key <- d[["key"]]
  url <- sprintf("%s/%s/%s.csv.gz", base, d[["db"]], d[["table"]])
  .odget(url, tmp, 600, sprintf("`%s`", key), fn)
  have <- .odhash(tmp, fn)
  if (have != want) {
    want <- .odsha(.odent(.odman(base, fn), key, base, fn))
  }
  if (is.na(want) || have != want) {
    .ckstop(
      fn,
      paste(
        "the file %s served for `%s` does not match the SHA-256 its",
        "manifest lists, also when read again (%s... against %s...): it",
        "was not used; the table may be being replaced, so try again in",
        "ten minutes"
      ),
      base,
      key,
      substr(have, 1L, 12L),
      substr(if (is.na(want)) "none" else want, 1L, 12L)
    )
  }
  if (!file.rename(tmp, path)) {
    .ckstop(fn, "the download could not be kept in %s", dir)
  }
  invisible(path)
}

# The kept file as a data frame, as served (column names unchanged); with
# cols, only those columns are read.
.odread <- function(path, cols, key, fn) {
  rd <- function(...) {
    tryCatch(
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                      encoding = "UTF-8", skipNul = TRUE, ...),
      error = function(e) {
        .ckstop(fn, "`%s` could not be read: %s", key, conditionMessage(e))
      }
    )
  }
  if (is.null(cols)) {
    return(rd())
  }
  have <- names(rd(nrows = 1L))
  i <- .odcol(cols, have, fn)
  cc <- rep("NULL", length(have))
  cc[i] <- NA
  rd(colClasses = cc)[match(i, which(is.na(cc)))]
}

# The positions of the columns cols among have: an exact name, else the
# one name equal to it without regard to case or underscores.
.odcol <- function(cols, have, fn) {
  if (!is.character(cols) || !length(cols) || anyNA(cols)) {
    .ckstop(fn, "`cols` must be column names, or NULL for all")
  }
  norm <- function(x) tolower(gsub("_", "", x, fixed = TRUE))
  vapply(cols, function(cl) {
    i <- match(cl, have)
    if (!is.na(i)) {
      return(i)
    }
    i <- which(norm(have) == norm(cl))
    if (length(i) > 1L) {
      .ckstop(
        fn, "`%s` matches more than one column: %s", cl, toString(have[i])
      )
    }
    if (!length(i)) {
      .ckstop(
        fn,
        "the table has no column `%s`; it has %s",
        cl,
        toString(have, width = 300L)
      )
    }
    i
  }, 1L, USE.NAMES = FALSE)
}

# url downloaded to dest, sent with a User-Agent the service accepts
# from a package, waiting up to wait seconds (or the session's longer
# timeout, restored afterwards). Any failure deletes dest and refuses
# in fn's name, saying what the service answered.
.odget <- function(url, dest, wait, what, fn) {
  to <- suppressWarnings(as.numeric(getOption("timeout", 60))[1L])
  old <- options(timeout = max(wait, to, na.rm = TRUE))
  on.exit(options(old), add = TRUE)
  ua <- c(
    "User-Agent" = paste0(
      "rmorie/1 (+https://crimlab.dev) crimlab/",
      utils::packageVersion("crimlab")
    )
  )
  rc <- tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = !interactive(),
                         headers = ua),
    error = function(e) e,
    warning = function(w) w
  )
  if (inherits(rc, "condition") || !identical(as.integer(rc), 0L)) {
    unlink(dest)
    m <- if (inherits(rc, "condition")) conditionMessage(rc) else ""
    .ckstop(fn, "%s", .odwhy(m, url, what, getOption("timeout")))
  }
  invisible(dest)
}

# A failed download in words, from download.file()'s message m.
.odwhy <- function(m, url, what, wait) {
  host <- sub("^(https?://[^/]*).*$", "\\1", url)
  st <- regmatches(m, regexec("status was '([0-9]{3})", m))[[1L]][2L]
  if (is.na(st)) {
    words <- c("Not Found" = "404", Forbidden = "403", Gone = "410",
               "Too Many" = "429", Unavailable = "503")
    hit <- vapply(names(words), grepl, NA, x = m, fixed = TRUE)
    st <- if (any(hit)) words[[which(hit)[1L]]] else NA_character_
  }
  if (grepl("timeout|timed out", m, ignore.case = TRUE)) {
    return(sprintf("%s did not answer within %s seconds; try again later",
                   host, format(wait)))
  }
  switch(
    if (is.na(st)) "" else st,
    "404" = sprintf("%s does not serve %s (404) at %s", host, what, url),
    "403" = sprintf("%s refused %s (403); update crimlab or try later",
                    host, what),
    "410" = sprintf(
      "%s does not serve %s (410): its publisher's terms do not allow it",
      host, what
    ),
    "429" = sprintf("%s has had too many requests; wait a minute", host),
    if (!is.na(st) && startsWith(st, "5")) {
      sprintf("%s is unavailable (%s); try again in a few minutes", host, st)
    } else {
      sprintf("%s could not be downloaded from %s (%s)", what, host,
              if (nzchar(m)) m else "no reason given")
    }
  )
}

# JSON text as lists (arrays too), with NULL for null. The strings stay as
# string literals; the text between them, once checked to hold only
# JSON's punctuation, numbers and literals, becomes list() calls, which
# are evaluated with nothing but list() and minus in scope.
# ponytail: \u escapes of characters outside the Basic Multilingual Plane
# (surrogate pairs) are refused with the manifest; decode them if a
# manifest ever carries one.
.odjson <- function(txt) {
  # bytes: character positions in a long UTF-8 text take minutes
  m <- gregexpr(r"("(?:[^"\\]++|\\.)*+")", txt, perl = TRUE, useBytes = TRUE)
  s <- regmatches(txt, m)[[1L]]
  g <- regmatches(txt, m, invert = TRUE)[[1L]]
  Encoding(s) <- "UTF-8"
  tok <- paste0(
    r"(^(?:\s|[\[\]{}:,]|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)",
    r"(|true|false|null)*$)"
  )
  esc <- r"(^"(?:[^\\]|\\["\\/bfnrtu])*"$)"
  ok <- all(grepl(tok, g, perl = TRUE, useBytes = TRUE)) &&
    all(grepl(esc, s, perl = TRUE, useBytes = TRUE))
  if (!ok) {
    stop("not JSON", call. = FALSE)
  }
  s <- gsub(r"((?<!\\)((?:\\\\)*)\\/)", "\\1/", s, perl = TRUE)
  g <- gsub(r"([\[{])", "list(", g, perl = TRUE)
  g <- gsub(r"([\]}])", ")", g, perl = TRUE)
  lit <- c(":" = "=", true = "TRUE", false = "FALSE", null = "NULL")
  for (k in names(lit)) g <- gsub(k, lit[[k]], g, fixed = TRUE)
  n <- length(s)
  e <- parse(
    text = paste(c(rbind(g[seq_len(n)], s), g[n + 1L]), collapse = ""),
    keep.source = FALSE,
    encoding = "UTF-8"
  )
  if (length(e) != 1L) {
    stop("not JSON", call. = FALSE)
  }
  eval(e[[1L]], list2env(list(list = list, "-" = `-`), parent = emptyenv()))
}
