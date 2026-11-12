# Time and peak memory of every exported function on one large input
# (cases in dev/bench.txt; n = 1e6 where n applies; scalar functions are timed over 1e4 calls,
# proof() over 1e3).
#   Rscript dev/bench.R LIB      prints: fn, median seconds of 5, peak MB
# Run it once with the reference build and once with the new one.
.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
library(crimlab)
grDevices::pdf(NULL)
set.seed(1)
n <- 1e6
u <- function(k, lo = 0.01, hi = 1) stats::runif(k, lo, hi)
x1 <- u(n, 1, 100)
x2 <- u(n, 1, 100)
p <- u(n, 0.01, 0.4)
g <- cbind(sample.int(30, n, TRUE), sample.int(30, n, TRUE))
h <- cbind(sample.int(30, n, TRUE), sample.int(30, n, TRUE))
fr <- matrix(0.02, n, 2)
gm <- 1000
cnt <- matrix(u(gm * gm, 1, 10), gm)
pg <- u(gm, 10, 100)
pm <- stats::setNames(rep(1e-6, 1e5), paste0("a", seq_len(1e5)))
# rep4(e, k): e k times, as a function so that the JIT compiles it
rep4 <- function(e, k = 1e4) {
  f <- eval(call("function", NULL, substitute(e)), parent.frame())
  for (i in seq_len(k)) f()
}
yb <- as.double(p > 0.2)
zl <- ifelse(x1 > 50, "a", "b")
ls1 <- sort(x1, decreasing = TRUE)
sh <- rep(1 / n, n)
ix2 <- 1 / x2
xs <- x1 + x2
s5 <- x1[1:5e5]
d5 <- ix2[1:5e5]
s4 <- x1[1:5e4]
d4 <- ix2[1:5e4]
tab <- utils::read.delim("dev/bench.txt", header = FALSE, quote = "")
cases <- stats::setNames(lapply(tab[[2]], function(s) {
  e <- str2lang(s)
  function() eval(e, globalenv())
}), tab[[1]])
only <- commandArgs(TRUE)[-1]
if (length(only)) cases <- cases[only]
for (nm in names(cases)) {
  f <- cases[[nm]]
  f()
  t <- vapply(1:5, function(i) system.time(f())[["elapsed"]], 0)
  invisible(gc(reset = TRUE))
  b <- sum(gc()[, 6])
  f()
  m <- sum(gc()[, 6]) - b
  cat(sprintf("%-8s %8.3f %8.1f\n", nm, stats::median(t), m))
}
