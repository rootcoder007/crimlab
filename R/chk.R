# Argument checks shared by every exported function. Each one stops with a
# message naming the argument, so the compiled core never sees a bad value.

.num <- function(x, nm, len = NULL, lo = -Inf, open = FALSE) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("`", nm, "` must be finite numbers with no missing values.",
         call. = FALSE)
  }
  if (!is.null(len) && length(x) != len) {
    stop("`", nm, "` must have length ", len, ", not ", length(x), ".",
         call. = FALSE)
  }
  bad <- if (open) x <= lo else x < lo
  if (any(bad)) {
    stop("`", nm, "` must be ", if (open) "greater than " else "at least ",
         lo, ".", call. = FALSE)
  }
  as.double(x)
}

.one <- function(x, nm, lo, hi) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= lo || x >= hi) {
    stop("`", nm, "` must be one number strictly between ", lo, " and ", hi,
         ".", call. = FALSE)
  }
  as.double(x)
}

# The position of a group given by name or by number, one-based.
.ref <- function(ref, nms, n) {
  if (is.character(ref) && length(ref) == 1L && !is.na(ref)) {
    i <- match(ref, nms)
    if (is.na(i)) stop("`ref` names no group: ", ref, ".", call. = FALSE)
    return(i)
  }
  if (is.numeric(ref) && length(ref) == 1L && !is.na(ref) &&
        ref == round(ref) && ref >= 1 && ref <= n) {
    return(as.integer(ref))
  }
  stop("`ref` must be one group name or a position from 1 to ", n, ".",
       call. = FALSE)
}

.grp <- function(x) {
  if (is.null(names(x))) as.character(seq_along(x)) else names(x)
}
