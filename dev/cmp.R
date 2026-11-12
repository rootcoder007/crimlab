# Old-versus-new comparison: each call in dev/cmp.txt (one per line; r()
# draws seeded uniforms, big holds 1000 ratios) is run against an installed
# build; values, attributes in order, messages, classes and calls compared.
#   Rscript dev/cmp.R run LIB OUT.rds
#   Rscript dev/cmp.R diff OLD.rds NEW.rds
a <- commandArgs(TRUE)
one <- function(s) {
  w <- list()
  v <- withCallingHandlers(
    tryCatch(eval(str2lang(s), globalenv()), error = function(e) {
      structure(list(class(e), conditionMessage(e), conditionCall(e)),
                class = "cmperr")
    }),
    warning = function(x) {
      w[[length(w) + 1L]] <<- list(conditionMessage(x), conditionCall(x))
      invokeRestart("muffleWarning")
    }
  )
  list(v, w)
}
near <- function(x, y) {
  if (is.list(x) && is.list(y) && length(x) == length(y) &&
        identical(attributes(x), attributes(y))) {
    return(all(mapply(near, x, y)))
  }
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
        !identical(attributes(x), attributes(y)) ||
        !identical(is.na(x), is.na(y))) {
    return(identical(x, y))
  }
  k <- !is.na(x)
  all(x[k] == y[k] | abs(x[k] - y[k]) <= 1e-12 * pmax(abs(x[k]), abs(y[k])))
}
if (a[1] == "run") {
  .libPaths(c(a[2], .libPaths()))
  library(crimlab)
  grDevices::pdf(NULL)
  set.seed(20261112)
  r <- function(n, lo = 0.01, hi = 1) stats::runif(n, lo, hi)
  big <- r(1000, 1e-3, 1e3)
  src <- readLines("dev/cmp.txt")
  saveRDS(list(src = src, out = lapply(src, one)), a[3])
} else {
  o <- readRDS(a[2])
  n <- readRDS(a[3])
  k <- c(0, 0, 0)
  for (i in seq_along(o$out)) {
    s <- if (identical(o$out[[i]], n$out[[i]], attrib.as.set = FALSE)) 1 else
      if (near(o$out[[i]], n$out[[i]])) 2 else 3
    k[s] <- k[s] + 1
    if (s > 1) cat(c("", "NEAR", "DIFF")[s], substr(o$src[i], 1, 100), "\n")
  }
  cat(sprintf("CMP identical %d, within 1e-12 %d, different %d\n", k[1],
              k[2], k[3]))
}
