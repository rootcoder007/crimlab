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

# The rows named by rows, permuted together three times, give the result
# of f(args): its numbers compared as a sorted set (cells and clusters come
# out in an order of their own, rows in the order given).
pinv <- function(f, args, rows) {
  key <- function(r) {
    v <- unname(rapply(list(r), as.numeric,
                       classes = c("numeric", "integer"), how = "unlist"))
    sort(signif(v[is.finite(v)], 10))
  }
  k0 <- key(do.call(f, args))
  n <- NROW(args[[rows[1L]]])
  set.seed(7)
  for (s in 1:3) {
    o <- sample(n)
    a <- args
    for (p in rows) {
      a[[p]] <- if (is.matrix(a[[p]])) a[[p]][o, , drop = FALSE] else a[[p]][o]
    }
    expect_equal(key(do.call(f, a)), k0, tolerance = 1e-10,
                 info = if (is.character(f)) f else "")
  }
}

# The cluster bootstrap's draws made all at once, g clusters by b draws
# from R's generator after set.seed(seed): what the package's block-wise
# .axdraws() must reproduce. The generator's state is put back.
.axdraw <- function(g, b, seed) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had) {
    assign(".Random.seed", old, envir = env)
  } else {
    rm(".Random.seed", envir = env)
  }, add = TRUE)
  set.seed(seed)
  .axdraw1(g, b)
}

