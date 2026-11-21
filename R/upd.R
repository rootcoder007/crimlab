# The notice of a newer release, when the package is attached.

#' Ask whether a newer release exists
#'
#' Asks the project's release page once, now, for the number of the
#' latest release, and says whether it is newer than the installed
#' package. Nothing is sent but that request (a limit of two seconds), and
#' nothing is kept.
#'
#' The package never asks on its own: no request is made when it is
#' attached unless a notice has been turned on (see the section "Notice
#' of a newer release" of [crimlab-package]).
#'
#' @return Invisibly, a list: `current` (the installed version), `latest`
#'   (the latest release, `NA` when the page could not be reached) and
#'   `newer` (`TRUE` when `latest` is above `current`). One line is
#'   printed.
#' @examples
#' \dontrun{
#' # needs the internet
#' updchk()
#' }
#' @export
updchk <- function() {
  cur <- as.character(utils::packageVersion("crimlab"))
  v <- updget()
  ok <- !is.na(v) && grepl("^[0-9]+(\\.[0-9]+)*$", v)
  new <- ok && utils::compareVersion(v, cur) > 0L
  cat(
    if (!ok) {
      "crimlab: the release page could not be reached\n"
    } else if (new) {
      sprintf("crimlab %s is available (this is %s): %s\n", v, cur, updurl)
    } else {
      sprintf("crimlab %s is the latest release\n", cur)
    }
  )
  invisible(list(
    current = cur,
    latest = if (ok) v else NA_character_,
    newer = new
  ))
}

# The address that answers with the address of the latest release.
updurl <- "https://github.com/rootcoder007/crimlab/releases/latest"

# The version in the header lines of that answer, or NA.
updver <- function(h) {
  l <- grep("^location:", h, ignore.case = TRUE, value = TRUE)
  m <- regmatches(l, regexec("/releases/tag/v([0-9]+(\\.[0-9]+)*)\\s*$", l))
  v <- unlist(lapply(m, function(x) x[2L]))
  v <- v[!is.na(v)]
  if (length(v) == 0L) NA_character_ else v[[1L]]
}

# The version of the latest release, or NA when it cannot be had: two
# seconds at most, and no error or warning leaves this function.
updget <- function(url = updurl) {
  op <- options(timeout = 2)
  on.exit(options(op), add = TRUE)
  tryCatch(
    suppressWarnings(updver(curlGetHeaders(url, redirect = FALSE))),
    error = function(e) NA_character_
  )
}

# Is the notice turned off. It is, unless it has been asked for:
# options(crimlab.update = TRUE), or the variable CRIMLAB_UPDATE_CHECK set to
# 1, true, yes or on. A session that is not interactive never asks, and
# options(crimlab.update = FALSE) or the variable CRIMLAB_NO_UPDATE_CHECK
# turn it off whatever else is set.
updoff <- function(live = interactive()) {
  on <- isTRUE(getOption("crimlab.update")) ||
    tolower(Sys.getenv("CRIMLAB_UPDATE_CHECK")) %in%
      c("1", "true", "yes", "on")
  !live || !on || nzchar(Sys.getenv("CRIMLAB_NO_UPDATE_CHECK")) ||
    isFALSE(getOption("crimlab.update"))
}

# The text of the notice, or NULL. The answer and its time are kept in
# `file`; the question is asked again after a day, by `get` (a test
# gives its own) or by updget().
updmsg <- function(cur, file, now = as.numeric(Sys.time()), get = NULL) {
  old <- tryCatch(
    scan(file, what = "", n = 2L, quiet = TRUE),
    error = function(e) character(),
    warning = function(w) character()
  )
  t <- suppressWarnings(as.numeric(old[1L]))
  v <- old[2L]
  if (is.na(t) || now < t || now - t >= 86400) {
    v <- if (is.null(get)) updget() else get()
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      writeLines(paste(format(now, scientific = FALSE), v), file),
      error = function(e) NULL,
      warning = function(w) NULL
    )
  }
  ok <- !is.na(v) && grepl("^[0-9]+(\\.[0-9]+)*$", v)
  if (!ok || utils::compareVersion(v, cur) <= 0L) {
    return(NULL)
  }
  paste0(
    "crimlab ", v, " is available (this is ", cur, "): ", updurl,
    "\noptions(crimlab.update = FALSE) stops this notice."
  )
}

# The notice for this installation, or NULL.
updrun <- function() {
  if (updoff()) {
    return(NULL)
  }
  updmsg(
    as.character(utils::packageVersion("crimlab")),
    file.path(tools::R_user_dir("crimlab", "cache"), "update")
  )
}

.onAttach <- function(libname, pkgname) {
  m <- updrun()
  if (!is.null(m)) {
    packageStartupMessage(m)
  }
}
