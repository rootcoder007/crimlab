# The tables of the project's open data service. The service's manifest,
# read at every call, says all there is to know about a table: its key,
# rows, columns, size and the SHA-256 of its gzip file as served. Nothing
# about any table is written here, so a publisher's new file needs no new
# release. A file is used only when its SHA-256 is the manifest's.

#' Tables of the open data service
#'
#' `odls()` lists the tables of the project's open data service and
#' `odtab()` loads one by its key, as the service serves it at the time of
#' the call. Everything about a table (its address, rows, columns and
#' SHA-256) is read from the service's manifest at run time, so a table
#' that its publisher has replaced is loaded as it is now, with no new
#' release of crimlab.
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
#' The file is read as CSV (RFC 4180): fields separated by commas, a field
#' in double quotes may hold commas, line breaks and doubled quotes, lines
#' end in LF or CRLF, a byte order mark is skipped. A blank line is a row
#' of a table of one column (its value missing), and is skipped in a wider
#' table, as `read.csv()` skips it. Bytes that
#' are not valid UTF-8 are kept as text, each written `<xx>` (its value in
#' hexadecimal), so no value is lost and none is in an invalid encoding.
#'
#' The columns are those of the file, named as served (an empty or a
#' repeated name too). `cols` picks some
#' of them, matched without regard to case or underscores (`"case_id"`
#' finds `"CaseID"`), so that a publisher's respelling of a name does not
#' break a call; an exact name is preferred.
#'
#' The values are those the file holds. Every column is read as text; an
#' empty field and the text `NA` are missing values. A column becomes
#' numbers only when every value in it is a plain decimal number that a
#' number holds exactly as written: digits with at most a leading minus,
#' a decimal point and an exponent, no leading zero, at most 15
#' significant digits, and not so large or so small that a number cannot
#' hold it (`"1e400"`, `"1e-400"` and `"1e-320"` stay text). So
#' identifiers such as `"007"` (not the same as `"7"`), long identifiers,
#' dates and codes such as `"F"` stay text, and no column becomes logical.
#' A column of whole numbers within R's integer range is integer, any other
#' numeric column double.
#'
#' Errors: a key the manifest does not list; a table listed without a
#' SHA-256; a file that does not match it; a line with more or fewer
#' fields than the header, a quote inside a field that is not in quotes,
#' a quoted field never closed, text after a closing quote or a NUL byte
#' (each naming the line); a table whose number of rows is not the one
#' the manifest lists; a column the table does not have, or a name that
#' matches more than one; a service that cannot be reached, answers with
#' an error or does not answer in time (60 seconds for the manifest, 600
#' for a table, or `getOption("timeout")` if longer).
#'
#' @param key A table's key, as `odls()` lists it.
#' @param cols The names of the columns to keep, or `NULL` (all).
#' @param refresh `TRUE` downloads the table again even when the cache
#'   holds a copy that matches.
#' @param base The address of the service, or of a mirror of it (a folder
#'   holding `manifest.json` and `<db>/<table>.csv.gz`, which may be a
#'   `file://` address). The default is `getOption("crimlab.data")`, or
#'   the project's service when that is not set.
#' @return `odls()` a data frame with a row per table: `key`, `rows`,
#'   `bytes` (of the gzip file), `sha256` (`NA` for a table that cannot be
#'   checked) and `columns`, a list of the column names. `odtab()` the
#'   table as a data frame.
#' @examples
#' \dontrun{
#' # needs the internet
#' tabs <- odls()
#' head(tabs[, c("key", "rows")])
#' x <- odtab(tabs$key[!is.na(tabs$sha256)][1L])
#' str(x)
#' }
#' @export
odtab <- function(
  key,
  cols = NULL,
  refresh = FALSE,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
) {
  fn <- "odtab"
  if (!is.character(key) || length(key) != 1L || is.na(key)) {
    .ckstop(fn, "`key` must be a single table key, as odls() lists them")
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
  .odread(path, cols, d, fn)
}

#' @rdname odtab
#' @export
odls <- function(
  base = getOption("crimlab.data", "https://data.crimlab.dev")
) {
  fn <- "odls"
  base <- .odbase(base, fn)
  ds <- Filter(is.list, .odman(base, fn)[["datasets"]])
  out <- data.frame(
    key = .odfld(ds, "key", NA_character_),
    rows = .odfld(ds, "rows", NA_real_),
    bytes = .odfld(ds, "bytes_gz", NA_real_),
    sha256 = vapply(ds, .odsha, ""),
    stringsAsFactors = FALSE
  )
  out$columns <- lapply(ds, function(d) as.character(unlist(d[["columns"]])))
  out
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
  .odget(sprintf("%s/manifest.json", base), f, 60, "the manifest", fn)
  m <- tryCatch(
    .odjson(readBin(f, "raw", file.size(f))),
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
  .ckstop(fn, "`%s` is not a table of %s: odls() lists them", key, base)
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

# One field of every entry, NA where it is missing or not a single value.
.odfld <- function(ds, f, na) {
  vapply(ds, function(d) {
    v <- d[[f]]
    if (length(v) != 1L || is.list(v)) {
      return(na)
    }
    suppressWarnings(as.vector(v, typeof(na)))
  }, na)
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

# The kept file as a data frame of the values it holds (column names
# unchanged), entry d of the manifest; with cols, only those columns.
.odread <- function(path, cols, d, fn) {
  key <- d[["key"]]
  r <- .Call(C_clab_r_csv, .odbytes(path, key, fn))
  .odbad(r[[1L]], key, fn)
  have <- r[[2L]]
  i <- if (is.null(cols)) seq_along(have) else .odcol(cols, have, fn)
  nr <- if (length(have)) length(r[[3L]][[1L]]) else 0L
  n <- d[["rows"]]
  if (is.numeric(n) && length(n) == 1L && is.finite(n) && n != nr) {
    .ckstop(fn, "`%s` has %d rows where its manifest lists %s: not used",
            key, nr, format(n))
  }
  # the names as served, empty and repeated ones too
  structure(lapply(r[[3L]][i], .odnum), names = have[i],
            class = "data.frame", row.names = .set_row_names(nr))
}

# The bytes of a (gzip) file, read in pieces of 16 MB.
.odbytes <- function(path, key, fn) {
  con <- tryCatch(gzfile(path, "rb"), error = function(e) {
    .ckstop(fn, "`%s` could not be read: %s", key, conditionMessage(e))
  })
  on.exit(close(con), add = TRUE)
  b <- list()
  repeat {
    x <- readBin(con, "raw", 16777216L)
    if (!length(x)) break
    b[[length(b) + 1L]] <- x
  }
  do.call(c, c(list(raw()), b))
}

# The refusal of a table that is not CSV, from the reader's status: the
# kind, the line where the record starts, its fields and the header's.
.odbad <- function(st, key, fn) {
  if (st[1L] == 0) {
    return(invisible())
  }
  if (st[1L] == 6) .ckstop(fn, "`%s` is empty: it has no header line", key)
  if (st[1L] == 1) {
    .ckstop(fn, "line %d of `%s` has %d fields where its header has %d",
            st[2L], key, st[3L], st[4L])
  }
  why <- c("has a quote inside a field that is not in quotes",
           "opens a quoted field that is never closed",
           "has text after the quote that closes a field",
           "holds a NUL byte")
  .ckstop(fn, paste("line %d of `%s`", why[st[1L] - 1]), st[2L], key)
}

# A column of text as numbers when every value is a plain decimal number
# a double holds exactly as written (see odtab()); else as it is.
.odnum <- function(v) {
  w <- v[!is.na(v)]
  if (!length(w) ||
        !all(grepl("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?$", w,
                   useBytes = TRUE))) {
    return(v)
  }
  n <- as.numeric(v)
  if (!.odexact(w, abs(n[!is.na(v)]))) {
    return(v)
  }
  if (all(grepl("^-?[0-9]+$", w)) && all(abs(n) <= .Machine$integer.max,
                                         na.rm = TRUE)) {
    return(as.integer(n))
  }
  n
}

# Whether the decimal numbers w are held exactly by their doubles (sizes
# nw): at most 15 significant digits, and a normal double unless written
# as 0; a value that underflows (1e-400 is 0, 1e-320 subnormal) or
# overflows is not.
.odexact <- function(w, nw) {
  sig <- nchar(sub("^0+", "", gsub("[-.]", "", sub("[eE].*$", "", w))))
  !any(sig > 15L) && all(is.finite(nw)) &&
    all(sig == 0L | nw >= .Machine$double.xmin)
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
      "crimlab/", utils::packageVersion("crimlab"),
      " (+https://crimlab.dev)"
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

# JSON text (or its bytes) as lists, by the reader of the core (src/json.h,
# strict RFC 8259, through src/bjsn.cpp): an array a list, an object a
# named list, a string UTF-8 text, a number a double, true and false
# logicals, null NULL. Nothing is evaluated; the reader is iterative, 10000
# levels deep, linear in the text.
.odjson <- function(txt) {
  if (is.character(txt)) txt <- charToRaw(txt)
  r <- .Call(C_clab_r_json, txt)
  if (!r[[1L]]) {
    stop("not JSON", call. = FALSE)
  }
  r[[2L]]
}
