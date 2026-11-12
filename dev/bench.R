# Time and peak memory of every exported function on one large input
# (n = 1e6 where n applies; scalar functions are timed over 1e4 calls,
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
rep4 <- function(f) for (i in seq_len(1e4)) f()
rep3 <- function(f) for (i in seq_len(1e3)) f()
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
cases <- list(
  vodtst = function() vodtst(x1, x2),
  orrisk = function() rep4(function() orrisk(3, 0.4, 0.5)),
  outtst = function() outtst(p, 0.2, w = x1),
  dspbm = function() dspbm(x1, x2, x1),
  dyadrt = function() dyadrt(cnt, pg),
  gtmp = function() gtmp(p, p, hom = FALSE, theta = 0.01),
  cllr = function() cllr(x1, ix2),
  lrpost = function() lrpost(x1, 0.01),
  mixlr = function() mixlr(pm, names(pm)[1:10], 2, names(pm)[1:5], 2),
  dbsrch = function() rep4(function() dbsrch(1e-6, 1e5, 3e7)),
  cllrmin = function() cllrmin(s5, d5),
  lrcal = function() lrcal(s5, d5),
  tipplt = function() tipplt(s5, d5),
  eceplt = function() eceplt(s4, d4),
  kinidx = function() kinidx(g, h, fr, "sibling"),
  patidx = function() patidx(g, h, fr, mother = g),
  clrsmp = function() clrsmp(x1, xs, x2, xs),
  expbnd = function() expbnd(x1, x2, 1.5),
  misrd = function() rep4(function() misrd(0.3, 0.1, 0.4, 0.9, 0.8)),
  rtmhs = function() rtmhs(x1, x2, 50, w = p),
  incap = function() incap(x1, 0.5, 2, share = sh),
  incapr = function() incapr(ls1, 10, 1000, 0.1),
  sellab = function() sellab(yb, x1 > 30, x2 > 30, w = p),
  losbnd = function() losbnd(x1, x2),
  sntbnd = function() sntbnd(yb, zl),
  hzdsel = function() rep4(function() hzdsel(0.3, 0.4, 0.1, c(0.9, 0.8),
                                             c(0.95, 0.9))),
  proof = function() rep3(function() proof(c("cllr", "gtmp")))
)
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
