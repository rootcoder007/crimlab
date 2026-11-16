# Argument checks shared by the search and rescue, air pollution,
# footprint and data description families: each refusal names the
# function, "fn: <message>.", as the command line gives it.

.ckstop <- function(fn, fmt, ...) {
  # a value quoted from the input may be text that is not valid UTF-8,
  # on which sprintf() fails: its bytes are written as <xx> instead
  a <- lapply(list(...), function(v) {
    if (is.character(v) && !all(validEnc(v))) {
      v <- iconv(v, "UTF-8", "UTF-8", sub = "byte")
    }
    v
  })
  stop(fn, ": ", do.call(sprintf, c(list(fmt), a)), ".", call. = FALSE)
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

# expr evaluated for the exported function fn: a refusal or warning from
# another exported function it calls ("ppmlfe: ...") is given again under
# fn's name, the one the user called.
.clas <- function(fn, expr) {
  ex <- getNamespaceExports("crimlab")
  re <- function(m) {
    k <- regexpr("^[a-z][a-z0-9]*: ", m)
    p <- if (k > 0L) substr(m, 1L, attr(k, "match.length") - 2L) else ""
    if (nzchar(p) && p != fn && p %in% ex) {
      paste0(fn, substring(m, nchar(p) + 1L))
    } else {
      m
    }
  }
  withCallingHandlers(
    expr,
    error = function(e) {
      m <- re(conditionMessage(e))
      if (!identical(m, conditionMessage(e))) stop(m, call. = FALSE)
    },
    warning = function(w) {
      m <- re(conditionMessage(w))
      if (!identical(m, conditionMessage(w))) {
        warning(m, call. = FALSE)
        invokeRestart("muffleWarning")
      }
    }
  )
}

# The choice argument `what` of the calling function, checked by .ckarg()
# against the choices its formals list (what match.arg(what) reads, with
# the package's refusal in place of match.arg()'s).
.ckmarg <- function(what, fn) {
  env <- parent.frame()
  ch <- eval(formals(sys.function(sys.parent()))[[what]], env)
  .ckarg(get(what, envir = env), ch, what, fn)
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

# A seed for R's own generator: one whole number from -2147483647 to
# 2^53 - 1; an integer when it fits.
.ckseed <- function(seed, what, fn) {
  if (
    !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || seed < -.Machine$integer.max || seed >= 2^53
  ) {
    .ckstop(
      fn,
      "`%s` must be a single whole number from -2147483647 to 2^53 - 1",
      what
    )
  }
  if (abs(seed) <= .Machine$integer.max) as.integer(seed) else as.numeric(seed)
}

# expr evaluated with R's generator set from seed (checked by .ckseed()),
# the caller's stream restored afterwards. Past 2147483647 the
# Mersenne-Twister state is that set.seed() makes from the seed's low 31
# bits, XORed with the state it makes from minus its high bits, shifted by
# one word, so that every seed up to 2^53 has a stream of its own.
.clseed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  if (is.integer(seed)) {
    set.seed(seed)
  } else {
    st <- function(s) {
      set.seed(s)
      get(".Random.seed", envir = globalenv())
    }
    lo <- st(as.integer(seed %% 2^31))
    hi <- st(-as.integer(seed %/% 2^31))
    w <- bitwXor(lo[3:626], hi[c(4:626, 3L)])
    w[is.na(w)] <- 1L
    assign(".Random.seed", c(lo[[1L]], 624L, w), envir = globalenv())
  }
  expr
}
