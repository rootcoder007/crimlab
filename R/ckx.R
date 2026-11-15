# Argument checks shared by the search and rescue, air pollution,
# footprint and data description families: each refusal names the
# function, "fn: <message>.", as the command line gives it.

.ckstop <- function(fn, fmt, ...) {
  stop(fn, ": ", sprintf(fmt, ...), ".", call. = FALSE)
}

# A numeric argument as a double vector (dimensions kept): numbers only,
# none missing or infinite, at least one; scalar: exactly one.
.cknum <- function(
  x,
  what,
  fn,
  scalar = TRUE,
  min = -Inf,
  max = Inf,
  positive = FALSE,
  whole = FALSE
) {
  if (is.null(x)) {
    .ckstop(fn, "`%s` must be numeric, got NULL", what)
  }
  if (is.logical(x) && !all(is.na(x))) {
    .ckstop(fn, "`%s` must be numeric, not logical", what)
  }
  if (
    is.factor(x) || is.list(x) || !(is.numeric(x) || is.logical(x))
  ) {
    .ckstop(fn, "`%s` must be numeric, got %s", what, class(x)[1L])
  }
  v <- as.numeric(x)
  if (!length(v)) {
    .ckstop(fn, "`%s` must have at least one value", what)
  }
  if (scalar && length(v) != 1L) {
    .ckstop(fn, "`%s` must be a single number", what)
  }
  if (anyNA(v) || !all(is.finite(v))) {
    .ckstop(fn, "`%s` must be finite, with no missing values", what)
  }
  if (positive && any(v <= 0)) {
    .ckstop(fn, "`%s` must be positive", what)
  }
  if (any(v < min | v > max)) {
    .ckstop(
      fn,
      "`%s` must be %s",
      what,
      if (is.finite(max)) {
        sprintf("between %s and %s", format(min), format(max))
      } else {
        sprintf("at least %s", format(min))
      }
    )
  }
  if (whole && any(v != round(v))) {
    .ckstop(fn, "`%s` must be a whole number", what)
  }
  if (!is.null(dim(x))) {
    dim(v) <- dim(x)
  }
  v
}

# One of the choices, matched partially; NULL or the whole default picks
# the first.
.ckarg <- function(arg, choices, what, fn) {
  if (is.null(arg) || identical(arg, choices)) {
    return(choices[[1L]])
  }
  ok <- is.character(arg) && length(arg) == 1L && !is.na(arg)
  i <- if (ok) pmatch(arg, choices, nomatch = 0L) else 0L
  if (i == 0L) {
    .ckstop(
      fn,
      "`%s` must be one of %s, not %s",
      what,
      toString(dQuote(choices, FALSE)),
      paste(utils::head(deparse(arg), 1L), collapse = "")
    )
  }
  choices[[i]]
}

# Text that is valid in its encoding (or anything that is not text).
.cktxt <- function(x, what, fn) {
  if (is.character(x) && !all(validEnc(x))) {
    .ckstop(
      fn,
      "`%s` holds text that is not valid %s",
      what,
      if (isTRUE(l10n_info()[["UTF-8"]])) "UTF-8" else "in this locale"
    )
  }
  x
}

# TRUE or FALSE.
.ckflag <- function(x, what, fn) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .ckstop(fn, "`%s` must be TRUE or FALSE", what)
  }
  x
}
