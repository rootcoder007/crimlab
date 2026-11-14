# Bad input, table-driven: each case is the message expected, then the
# arguments of one call to f.
#' @srrstats {G5.2b} errs() triggers each refusal and matches its message.
# Every message starts with the name of the function that raised it.
errs <- function(f, cases) {
  for (a in cases) {
    e <- expect_error(do.call(f, a[-1L]), a[[1L]], info = f)
    expect_match(conditionMessage(e), paste0("^", f, ": "), info = f)
  }
}

# Reference values from other packages, stored with the tests: the list
# written to ref/<name>.txt by dev/ref/<name>.R, whose first lines name
# the packages, versions and platform that computed it. With CLAB_REF set
# to a directory of regenerated fixtures (dev/ref/all.R; the ref
# workflow), a fixture found there is read instead.
refget <- function(name) {
  f <- paste0(name, ".txt")
  d <- Sys.getenv("CLAB_REF")
  if (nzchar(d) && file.exists(file.path(d, f))) {
    return(dget(file.path(d, f)))
  }
  dget(test_path("ref", f))
}

# The leaves of a stored list: its atomic vectors, each named by its path
# ("stat/gc", "rand/1u/btw"; unnamed elements by their position), lists
# and data frames walked into. dev/ref/cmp.R uses it too.
refleaf <- function(x, p = NULL) {
  if (!is.list(x)) {
    return(stats::setNames(list(x), paste(p, collapse = "/")))
  }
  k <- names(x)
  if (is.null(k)) k <- character(length(x))
  k[k == ""] <- which(k == "")
  unlist(lapply(seq_along(x), function(i) refleaf(x[[i]], c(p, k[i]))),
         recursive = FALSE)
}

# For each leaf path, the index of its row among the patterns pat (the
# path column of ref/tol.csv for one fixture): the longest pattern whose
# segments begin the path, a "*" segment matching any one; NA for none.
refrow <- function(pat, paths) {
  q <- strsplit(pat, "/", fixed = TRUE)
  vapply(strsplit(paths, "/", fixed = TRUE), function(s) {
    ok <- vapply(q, function(x) {
      length(x) <= length(s) && all(x == "*" | x == s[seq_along(x)])
    }, NA)
    if (!any(ok)) return(NA_integer_)
    i <- which(ok)
    i[which.max(lengths(q)[i])]
  }, 1L)
}
