# print, format and summary of the results of the scaling (class svfit)
# and count-dynamics (class tsfit) functions: one line per part. The print
# method of these and of every other result whose format method gives its
# lines (the Indigenous health and wellbeing results of R/igc.R, R/igf.R,
# R/igh.R and R/igm.R) is the one function .prfmt(), registered here for
# each class.

#' @export
format.svfit <- function(x, ...) {
  nm <- names(x)
  c(
    sprintf("%s:", class(x)[1L]),
    sprintf("  %-*s %s", max(nchar(nm)), nm, vapply(unclass(x), .svpart, ""))
  )
}

#' @export
summary.svfit <- function(object, ...) {
  nm <- names(object)
  num <- vapply(unclass(object), function(v) is.numeric(v) && length(v), NA)
  rng <- function(v, f) {
    if (is.numeric(v) && length(v) && !all(is.na(v))) f(v, na.rm = TRUE) else NA
  }
  s <- data.frame(
    part = nm,
    size = vapply(unclass(object), function(v) {
      if (is.null(dim(v))) length(v) else prod(dim(v))
    }, 1),
    min = vapply(unclass(object), rng, 1, f = min),
    max = vapply(unclass(object), rng, 1, f = max),
    stringsAsFactors = FALSE
  )
  s$numeric <- num
  structure(list(fn = class(object)[1L], parts = s),
            class = "summary.svfit")
}

#' @export
print.summary.svfit <- function(x, ...) {
  cat(x$fn, ": ", nrow(x$parts), " parts\n", sep = "")
  print(x$parts, row.names = FALSE)
  invisible(x)
}

#' @export
format.tsfit <- function(x, ...) format.svfit(x, ...)

#' @export
summary.tsfit <- function(object, ...) summary.svfit(object, ...)

.svpart <- function(v) {
  if (is.null(v)) {
    return("NULL")
  }
  if (is.list(v)) {
    return(sprintf("list of %d", length(v)))
  }
  d <- dim(v)
  if (!is.null(d)) {
    return(sprintf("%s array", paste(d, collapse = " x ")))
  }
  if (length(v) == 1L) {
    return(if (is.numeric(v)) format(signif(v, 5)) else format(v))
  }
  if (is.numeric(v) && length(v) <= 6L) {
    return(paste(format(signif(v, 4)), collapse = " "))
  }
  sprintf("%d %s", length(v), if (is.character(v)) "strings" else "values")
}

# print a result: the lines its format method gives, then the result,
# invisibly
.prfmt <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @export
print.cnsprof <- .prfmt

#' @export
print.profshr <- .prfmt

#' @export
print.cnscwb <- .prfmt

#' @export
print.trchlthgap <- .prfmt

#' @export
print.mhcont <- .prfmt

#' @export
print.cwbscr <- .prfmt

#' @export
print.connassoc <- .prfmt

#' @export
print.healchg <- .prfmt

#' @export
print.commidx <- .prfmt

#' @export
print.frmaud <- .prfmt

#' @export
print.sctidy <- .prfmt

#' @export
print.limfrom <- .prfmt

#' @export
print.indgap <- .prfmt

#' @export
print.exptrnd <- .prfmt

#' @export
print.trcprog <- .prfmt

#' @export
print.ordcmp <- .prfmt

#' @export
print.stdprev <- .prfmt

#' @export
print.welprof <- .prfmt

#' @export
print.holbal <- .prfmt

#' @export
print.dualcont <- .prfmt

#' @export
print.svfit <- .prfmt

#' @export
print.tsfit <- .prfmt
