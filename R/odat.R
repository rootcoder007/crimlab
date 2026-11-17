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
#' The columns are those of the file, named as served. `cols` picks some
#' of them, matched without regard to case or underscores (`"case_id"`
#' finds `"CaseID"`), so that a publisher's respelling of a name does not
#' break a call; an exact name is preferred.
#'
#' The values are those the file holds. Every column is read as text; an
#' empty field and the text `NA` are missing values. A column becomes
#' numbers only when every value in it is a plain decimal number that a
#' number holds exactly as written: digits with at most a leading minus,
#' a decimal point and an exponent, no leading zero, at most 15
#' significant digits. So identifiers such as `"007"` (not the same as
#' `"7"`), long identifiers, dates and codes such as `"F"` stay text, and
#' no column becomes logical. A column of whole numbers within R's
#' integer range is integer, any other numeric column double.
#'
#' Errors: a key the manifest does not list; a table listed without a
#' SHA-256; a file that does not match it; a line with more or fewer
#' fields than the header; a table whose number of rows is not the one
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
  rd <- function(...) {
    tryCatch(
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                      encoding = "UTF-8", skipNul = TRUE, fill = FALSE,
                      na.strings = c("", "NA"), ...),
      error = function(e) {
        .ckstop(fn, "`%s` could not be read: %s", key, conditionMessage(e))
      }
    )
  }
  .odrag(path, key, fn)
  have <- names(rd(nrows = 1L, colClasses = "character"))
  i <- if (is.null(cols)) seq_along(have) else .odcol(cols, have, fn)
  cc <- rep("NULL", length(have))
  cc[i] <- "character"
  x <- rd(colClasses = cc)[match(i, which(cc != "NULL"))]
  n <- d[["rows"]]
  if (is.numeric(n) && length(n) == 1L && is.finite(n) && n != nrow(x)) {
    .ckstop(fn, "`%s` has %d rows where its manifest lists %s: not used",
            key, nrow(x), format(n))
  }
  x[] <- lapply(x, .odnum)
  x
}

# A file whose lines do not all have as many fields as its header (a
# field in quotes may span lines) is refused, naming the first such line:
# read.csv() would wrap a long line into a row of its own.
.odrag <- function(path, key, fn) {
  nf <- utils::count.fields(path, sep = ",", quote = "\"", comment.char = "",
                            blank.lines.skip = FALSE)
  bad <- which(!is.na(nf) & nf != 0L & nf != nf[1L])
  if (length(bad)) {
    .ckstop(fn, "line %d of `%s` has %d fields where its header has %d",
            bad[1L], key, nf[bad[1L]], nf[1L])
  }
}

# A column of text as numbers when every value is a plain decimal number
# a double holds exactly as written (see odtab()); else as it is.
.odnum <- function(v) {
  w <- v[!is.na(v)]
  if (!length(w) ||
        !all(grepl("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?$", w))) {
    return(v)
  }
  sig <- nchar(sub("^0+", "", gsub("[-.]", "", sub("[eE].*$", "", w))))
  n <- as.numeric(v)
  if (any(sig > 15L) || !all(is.finite(n[!is.na(v)]))) {
    return(v)
  }
  if (all(grepl("^-?[0-9]+$", w)) && all(abs(n) <= .Machine$integer.max,
                                         na.rm = TRUE)) {
    return(as.integer(n))
  }
  n
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

# JSON text as lists (arrays too), with NULL for null. The strings stay as
# string literals; the text between them, once checked to hold only
# JSON's punctuation, numbers and literals, becomes calls of .odarr() and
# .odobj(), a key naming its value's argument ("key" = value; the empty
# key, which R cannot name an argument, as "" := value), evaluated with
# nothing but those, `:=` and minus in scope. No pattern repeats per
# character of a string, so a string of any length is read.
.odjson <- function(txt) {
  # bytes: character positions in a long UTF-8 text take minutes
  m <- gregexpr(r"("[^"\\]*+(?:\\.[^"\\]*+)*+")", txt, perl = TRUE,
                useBytes = TRUE)
  s <- regmatches(txt, m)[[1L]]
  g <- regmatches(txt, m, invert = TRUE)[[1L]]
  Encoding(s) <- "UTF-8"
  # between the strings only JSON's punctuation, numbers and literals,
  # and in them every backslash begins one of its escapes: nothing is
  # left once they are taken out
  tok <- r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null|[\s\[\]{}:,]+)"
  esc <- r"(\\(?:["\\/bfnrt]|u[0-9A-Fa-f]{4}))"
  ok <- !any(nzchar(gsub(tok, "", g, perl = TRUE, useBytes = TRUE))) &&
    !any(grepl("\\", gsub(esc, "", s, perl = TRUE, useBytes = TRUE),
               fixed = TRUE))
  if (!ok) {
    stop("not JSON", call. = FALSE)
  }
  s <- .odesc(s)
  g <- gsub("[", ".odarr(", g, fixed = TRUE)
  g <- gsub("{", ".odobj(", g, fixed = TRUE)
  g <- gsub("[]}]", ")", g, perl = TRUE)
  n <- length(s)
  lit <- c(true = "TRUE", false = "FALSE", null = "NULL")
  for (k in names(lit)) g <- gsub(k, lit[[k]], g, fixed = TRUE)
  # the gap after an empty key starts with its colon
  e0 <- c(FALSE, s == r"("")")
  g[e0] <- sub(":", ":=", g[e0], fixed = TRUE)
  g[!e0] <- gsub(":", "=", g[!e0], fixed = TRUE)
  e <- parse(
    text = paste(c(rbind(g[seq_len(n)], s), g[n + 1L]), collapse = ""),
    keep.source = FALSE,
    encoding = "UTF-8"
  )
  if (length(e) != 1L) {
    stop("not JSON", call. = FALSE)
  }
  # with no empty key, no pair can arise: the quicker checks suffice
  f <- if (any(e0)) list(.odarr0, .odobj0) else list(.odarr, .odobj)
  r <- eval(e[[1L]], list2env(
    list(.odarr = f[[1L]], .odobj = f[[2L]], ":=" = .odkv, "-" = `-`),
    parent = emptyenv()
  ))
  if (.odis(list(r))) {
    stop("not JSON", call. = FALSE)
  }
  r
}

# The empty key and its value, as a pairlist (which no JSON value is); an
# array; an object, its keys as names. Each refuses what JSON's grammar
# does not allow: a key in an array, a value without a key in an object.
# .odarr0() and .odobj0() take the pairs of empty keys too.
.odkv <- function(k, v) {
  if (!identical(k, "") || (is.pairlist(v) && !is.null(v))) {
    stop("not JSON", call. = FALSE)
  }
  pairlist(k, v)
}

.odarr <- function(...) {
  a <- list(...)
  if (any(nzchar(names(a)))) {
    stop("not JSON", call. = FALSE)
  }
  a
}

.odobj <- function(...) {
  a <- list(...)
  nm <- names(a)
  if (length(a) && (is.null(nm) || !all(nzchar(nm)))) {
    stop("not JSON", call. = FALSE)
  }
  a
}

.odarr0 <- function(...) {
  a <- .odarr(...)
  if (any(.odis(a))) {
    stop("not JSON", call. = FALSE)
  }
  a
}

.odobj0 <- function(...) {
  a <- list(...)
  kv <- .odis(a)
  nm <- names(a)
  if (is.null(nm)) {
    nm <- character(length(a))
  }
  if (any(kv == nzchar(nm))) {
    stop("not JSON", call. = FALSE)
  }
  a[kv] <- lapply(a[kv], .subset2, 2L)
  if (length(a)) {
    names(a) <- nm
  }
  a
}

# Which elements of the list a are pairs.
.odis <- function(a) vapply(a, is.pairlist, NA) & !vapply(a, is.null, NA)

# JSON string literals as R's: \/ is /, and the \u escapes R does not
# take (a pair for a character past U+FFFF, a lone half of one, U+0000)
# are written as R's \U{} or as U+FFFD.
.odesc <- function(s) {
  pre <- r"((?<!\\)((?:\\\\)*))"
  s <- gsub(paste0(pre, r"(\\/)"), "\\1/", s, perl = TRUE, useBytes = TRUE)
  sur <- paste0(pre, r"(\\u[dD][89abAB][0-9a-fA-F]{2})",
                r"(\\u[dD][c-fC-F][0-9a-fA-F]{2})")
  hit <- grepl(sur, s, perl = TRUE, useBytes = TRUE)
  if (any(hit)) {
    m <- gregexpr(sur, s[hit], perl = TRUE, useBytes = TRUE)
    regmatches(s[hit], m) <- lapply(regmatches(s[hit], m), function(x) {
      h <- strtoi(substr(x, nchar(x) - 9L, nchar(x) - 6L), 16L)
      l <- strtoi(substr(x, nchar(x) - 3L, nchar(x)), 16L)
      cp <- 65536L + (h - 55296L) * 1024L + (l - 56320L)
      paste0(substr(x, 1L, nchar(x) - 12L), sprintf("\\U{%X}", cp))
    })
  }
  gsub(paste0(pre, r"(\\u(?:[dD][89a-fA-F][0-9a-fA-F]{2}|0000))"),
       "\\1\\\\ufffd", s, perl = TRUE, useBytes = TRUE)
}
