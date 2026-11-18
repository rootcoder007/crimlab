# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# For each help page: runs its example and names the object (if any) that the
# package's own plot() method draws without error, and how long the example took.
#   Rscript vigplot.R <package> <rd.json> <out.tsv>
a <- commandArgs(TRUE)
pkg <- a[1L]
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
ns <- asNamespace(pkg)
d <- jsonlite::fromJSON(a[2L], simplifyVector = FALSE)
grDevices::pdf(NULL)
out <- character(0)
for (p in d) {
  e <- new.env(parent = globalenv())
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    set.seed(1)
    suppressWarnings(suppressMessages(utils::capture.output(
      eval(parse(text = p$ex), e))))
    TRUE
  }, error = function(err) conditionMessage(err))
  el <- proc.time()[["elapsed"]] - t0
  v <- ""
  if (isTRUE(ok)) {
    al <- unlist(p$aliases)
    cand <- character(0)
    for (n in ls(e)) {
      x <- get(n, e)
      cl <- class(x)
      m <- NULL
      for (k in cl) {
        m <- utils::getS3method("plot", k, optional = TRUE, envir = ns)
        if (!is.null(m) && identical(environment(m), ns)) break
        m <- NULL
      }
      if (is.null(m)) next
      good <- tryCatch({
        suppressWarnings(plot(x))
        TRUE
      }, error = function(err) FALSE)
      if (good) cand <- c(cand, n)
    }
    if (length(cand)) {
      # the object the example makes last, as far as the text tells
      pos <- vapply(cand, function(n) {
        g <- gregexpr(paste0("(^|\\n)\\s*", n, "\\s*<-"), p$ex)[[1L]]
        if (g[1L] < 0) 0L else max(g)
      }, 0L)
      v <- cand[which.max(pos)]
    }
  }
  out[length(out) + 1L] <- paste(p$file, if (isTRUE(ok)) "ok" else "ERR", v,
                                 sprintf("%.2f", el),
                                 if (isTRUE(ok)) "" else gsub("[\t\n]", " ", ok),
                                 sep = "\t")
}
writeLines(out, a[3L])
cat(length(out), "pages;", sum(grepl("\tERR\t", out)), "errors;",
    sum(vapply(strsplit(out, "\t"), function(x) nzchar(x[3L]), NA)), "with a plot\n")
