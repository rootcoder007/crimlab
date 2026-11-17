# Argument checks shared by the search and rescue, air pollution,
# footprint and data description families: each refusal names the
# function, "fn: <message>.", as the command line gives it.

.ckstop <- function(fn, fmt, ...) {
  # a value quoted from the input may be text that is not valid UTF-8,
  # on which sprintf() fails: its bytes are written as <xx> instead
  a <- lapply(list(...), .u8)
  stop(fn, ": ", do.call(sprintf, c(list(fmt), a)), ".", call. = FALSE)
}

# Text with every byte that is not valid in its encoding written as <xx>
# (anything that is not text as it is).
.u8 <- function(x) {
  if (is.character(x) && !all(validEnc(x))) {
    x <- iconv(x, "UTF-8", "UTF-8", sub = "byte")
  }
  x
}

# stop(), warning() and message() of the package: base R's, with the
# message made valid text by .u8(). A message that quotes the input may
# hold its invalid bytes, and substr(), grepl() or nchar() on such a
# condition message fail; defined once here, no refusal anywhere in the
# package can carry them. A condition object passes with its message
# made valid the same way.
stop <- function(..., call. = TRUE, domain = NULL) {
  cl <- sys.call(-1L)
  base::stop(.u8cnd(list(...), simpleError, call., domain, cl))
}

warning <- function(..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
                    domain = NULL) {
  cl <- sys.call(-1L)
  cnd <- .u8cnd(list(...), simpleWarning, call., domain, cl)
  # base R prints a warning at once with immediate.; a condition is printed
  # at once when the option warn is 1 (noBreaks. only changes how base R
  # wraps that line, and is taken for the same arguments)
  if (isTRUE(immediate.) && getOption("warn") < 1L) {
    op <- options(warn = 1L)
    on.exit(options(op), add = TRUE)
  }
  base::warning(cnd)
}

message <- function(..., domain = NULL, appendLF = TRUE) {
  a <- list(...)
  if (length(a) == 1L && inherits(a[[1L]], "condition")) {
    a[[1L]][["message"]] <- .u8(conditionMessage(a[[1L]]))
    return(base::message(a[[1L]]))
  }
  # the call is this one, message(...), as base R's
  m <- .makeMessage(..., domain = domain, appendLF = appendLF)
  base::message(simpleMessage(.u8(m), sys.call()))
}

# The condition stop() or warning() signals: the one given, or one built
# by mk from the arguments, with the call when call. is TRUE.
.u8cnd <- function(a, mk, call., domain, cl) {
  if (length(a) == 1L && inherits(a[[1L]], "condition")) {
    a[[1L]][["message"]] <- .u8(conditionMessage(a[[1L]]))
    return(a[[1L]])
  }
  m <- do.call(.makeMessage, c(a, list(domain = domain)))
  mk(.u8(m), if (call.) cl)
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
  # a scalar is one plain number: a 1 x 1 matrix keeps no dimensions,
  # which R's arithmetic would recycle with a warning
  if (!scalar && !is.null(dim(x))) {
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
  # by bytes: the message may quote text that is not valid in its encoding
  re <- function(m) {
    k <- regexpr("^[a-z][a-z0-9]*: ", m, useBytes = TRUE)
    if (k < 0L) {
      return(m)
    }
    b <- charToRaw(m)
    n <- attr(k, "match.length") - 2L
    p <- rawToChar(b[seq_len(n)])
    if (p != fn && p %in% ex) paste0(fn, rawToChar(b[-seq_len(n)])) else m
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

# Every text column of a data frame (a factor's labels, and the names
# too) valid text, as .cktxt() asks of one argument.
.cktxtdf <- function(d, what, fn) {
  v <- c(list(names(d)),
         lapply(d, function(x) if (is.factor(x)) levels(x) else x))
  for (x in v) {
    if (is.character(x)) .cktxt(x, what, fn)
  }
  d
}

# One file path (NULL passes): valid text, and no longer than a path can
# be (4096 bytes, PATH_MAX on Linux; R's path.expand() warns past it).
.ckpath <- function(file, fn, what = "file") {
  .cktxt(file, what, fn)
  if (!is.null(file) && nchar(file, "bytes") > 4096L) {
    .ckstop(fn, "`%s` is longer than a file path can be (4096 bytes)", what)
  }
  file
}

# A file to write (NULL passes): a path as .ckpath() takes it, not empty,
# in a folder that exists (where R's file() would warn, then fail in its
# own words).
.ckout <- function(file, fn, what = "file") {
  .ckpath(file, fn, what)
  if (is.null(file)) {
    return(file)
  }
  if (!nzchar(file)) .ckstop(fn, "`%s` must be a file path, not empty", what)
  if (!dir.exists(dirname(path.expand(file)))) {
    .ckstop(fn, "`%s` is in a folder that does not exist: %s", what,
            dirname(file))
  }
  file
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

# One text key per row of the vectors given (of one length, or 1): two
# rows share a key exactly when they agree in every vector. Each value is
# written with its backslashes doubled and its carriage returns as \r, a
# missing value as \N, and the values are joined by a carriage return, so
# no label can imitate the join (as with paste(sep = "\r"), where "a\rb"
# and "a", "b" met).
.clkey <- function(...) {
  v <- lapply(list(...), function(x) {
    s <- gsub("\\", "\\\\", as.character(x), fixed = TRUE)
    s <- gsub("\r", "\\r", s, fixed = TRUE)
    s[is.na(x)] <- "\\N"
    s
  })
  do.call(paste, c(unname(v), sep = "\r"))
}
