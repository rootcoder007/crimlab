# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# For each help page: runs its example and names what the package's own
# plot() method draws without error (an object the example saves, or else
# the last call whose result is only printed), and how long the example
# took.
#   Rscript vigplot.R <package> <rd.json> <out.tsv>
# out.tsv: file, ok or ERR, the object to plot, seconds, the error, and the
# text of the call to plot when no object was saved (newlines as \n).
a <- commandArgs(TRUE)
pkg <- a[1L]
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
ns <- asNamespace(pkg)
d <- jsonlite::fromJSON(a[2L], simplifyVector = FALSE)
grDevices::pdf(NULL)

# does the package draw x with a method of its own
draws <- function(x) {
  m <- NULL
  for (k in class(x)) {
    m <- utils::getS3method("plot", k, optional = TRUE, envir = ns)
    if (!is.null(m) && identical(environment(m), ns)) break
    m <- NULL
  }
  if (is.null(m)) return(FALSE)
  tryCatch({
    suppressWarnings(plot(x))
    TRUE
  }, error = function(err) FALSE)
}

out <- character(0)
for (p in d) {
  e <- new.env(parent = globalenv())
  t0 <- proc.time()[["elapsed"]]
  call <- ""
  ok <- tryCatch({
    set.seed(1)
    ex <- parse(text = p$ex, keep.source = TRUE)
    sr <- attr(ex, "srcref")
    for (i in seq_along(ex)) {
      v <- NULL
      suppressWarnings(suppressMessages(utils::capture.output(
        v <- withVisible(eval(ex[[i]], e)))))
      saved <- is.call(ex[[i]]) &&
        as.character(ex[[i]][[1L]])[1L] %in% c("<-", "=", "<<-")
      # a call that only prints its result: the last one that can be drawn
      if (!saved && v$visible && is.call(ex[[i]]) && draws(v$value)) {
        call <- paste(as.character(sr[[i]]), collapse = "\n")
      }
    }
    TRUE
  }, error = function(err) conditionMessage(err))
  el <- proc.time()[["elapsed"]] - t0
  v <- ""
  if (isTRUE(ok)) {
    cand <- character(0)
    for (n in ls(e)) if (draws(get(n, e))) cand <- c(cand, n)
    if (length(cand)) {
      # the object the example makes last, as far as the text tells
      pos <- vapply(cand, function(n) {
        g <- gregexpr(paste0("(^|\\n)\\s*", n, "\\s*<-"), p$ex)[[1L]]
        if (g[1L] < 0) 0L else max(g)
      }, 0L)
      v <- cand[which.max(pos)]
      call <- ""
    }
  }
  out[length(out) + 1L] <- paste(p$file, if (isTRUE(ok)) "ok" else "ERR", v,
                                 sprintf("%.2f", el),
                                 if (isTRUE(ok)) "" else gsub("[\t\n]", " ", ok),
                                 gsub("\n", "\\n", call, fixed = TRUE),
                                 sep = "\t")
}
writeLines(out, a[3L])
f <- strsplit(out, "\t")
cat(length(out), "pages;", sum(grepl("\tERR\t", out)), "errors;",
    sum(vapply(f, function(x) nzchar(x[3L]), NA)), "with a saved object to plot,",
    sum(vapply(f, function(x) length(x) > 5L && nzchar(x[6L]), NA)),
    "with a call to plot\n")
