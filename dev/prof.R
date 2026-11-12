# One call of crimlab, repeated, for a profiler: the input is made first,
# then the call runs k times.
#   valgrind --tool=callgrind Rscript dev/prof.R LIB CASE [k]
# CASE: vodtst, dspbm, incap, outtst, sntbnd, kinidx, ntdeg, small (each a
# 1e5-group, 1e6-value or 1e5-node input, named; small: the same calls at
# 5 groups, unnamed).
a <- commandArgs(TRUE)
.libPaths(c(a[1], .libPaths()))
library(crimlab)
k <- if (length(a) > 2) as.integer(a[3]) else 10L
set.seed(1)
n <- 1e5
g <- paste0("g", seq_len(n))
v <- function() stats::setNames(stats::runif(n, 10, 1e3), g)
f <- switch(a[2],
  vodtst = { d <- v(); e <- v(); function() vodtst(d, e) },
  dspbm = { p <- v(); cn <- v(); fo <- v(); function() dspbm(p, cn, fo) },
  incap = { l <- v(); sh <- stats::runif(n); sh <- sh / sum(sh)
            function() incap(unname(l), 0.1, 2, share = sh) },
  outtst = { q <- stats::runif(1e6); w <- stats::runif(1e6)
             names(q) <- names(w) <- paste0("p", seq_len(1e6))
             function() outtst(q, 0.5, w) },
  sntbnd = { z <- sample(c("a", "b"), 1e6, TRUE); y <- stats::rbinom(1e6, 1, 0.4)
             names(y) <- names(z) <- paste0("p", seq_len(1e6))
             function() sntbnd(y, z, "b") },
  kinidx = { s <- function() t(replicate(n, sample(6, 2, TRUE)))
             g1 <- s(); g2 <- s(); fr <- matrix(stats::runif(2 * n, 0.05, 0.3), n)
             rownames(g1) <- rownames(g2) <- rownames(fr) <- g
             function() kinidx(g1, g2, fr, "sibling") },
  ntdeg = { fr <- sample(g, 3 * n, TRUE); to <- sample(g, 3 * n, TRUE)
            function() ntdeg(fr, to, nodes = g) },
  small = { x <- stats::runif(5); y <- stats::runif(5)
            function() { vodtst(x, y); dspbm(x, y, x); clrsmp(x, y, x, y)
                         incap(x, 0.1, 2); misrd(0.3, 0.1, 0.4, 0.8, 0.9) } })
invisible(f())
for (i in seq_len(k)) invisible(f())
