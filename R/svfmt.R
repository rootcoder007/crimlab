# print, format and summary of the results of the scaling (class svfit)
# functions: one line per part.

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

#' @export
format.svfit <- function(x, ...) {
  nm <- names(x)
  c(
    sprintf("%s:", class(x)[1L]),
    sprintf("  %-*s %s", max(nchar(nm)), nm, vapply(unclass(x), .svpart, ""))
  )
}

#' @export
print.svfit <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @export
summary.svfit <- function(object, ...) {
  nm <- names(object)
  num <- vapply(unclass(object), function(v) is.numeric(v) && length(v), NA)
  rng <- function(v, f) {
    if (is.numeric(v) && length(v) && any(!is.na(v))) f(v, na.rm = TRUE) else NA
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

