# Reference values and checks for src/sens.cpp (clab_rswp, clab_rshl).
# Field package: rbounds 2.2 (CRAN), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/rbounds.R
library(rbounds)
cat("rbounds", format(packageVersion("rbounds")), "\n")
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")

# deterministic pair differences, rebuilt in src/t_sens.cpp: multiples of
# 1/16 (exact in binary), with ties and two zero differences
i <- 0:59
d <- (((i * 37) %% 101) - 35) / 16
cat("zero differences:", sum(d == 0), " tied |d| values:", sum(duplicated(abs(d[d != 0]))), "\n")
g <- seq(1, 3, by = 0.5)
ps <- psens(d, rep(0, 60), Gamma = 3, GammaInc = 0.5)$bounds
cat("\n// psens(d, 0, Gamma = 3, GammaInc = 0.5)$bounds (rounded to 4 digits by rbounds)\n")
cat(sprintf("const double ref_plo[5] = {%s};\n", f17(ps[, 2])))
cat(sprintf("const double ref_pup[5] = {%s};\n", f17(ps[, 3])))
# the same normal-approximation bounds without the rounding, from the
# formulas in rbounds' psens (base R)
dd <- d[d != 0]; r <- rank(abs(dd)); Tp <- sum(r[dd > 0])
pu <- sapply(g, function(G) { p <- G / (1 + G); pnorm((Tp - p * sum(r)) / sqrt(p * (1 - p) * sum(r^2)), lower.tail = FALSE) })
pl <- sapply(g, function(G) { p <- 1 / (1 + G); pnorm((Tp - p * sum(r)) / sqrt(p * (1 - p) * sum(r^2)), lower.tail = FALSE) })
cat("// unrounded, base R pnorm(lower.tail = FALSE) on the same formulas\n")
cat(sprintf("const double ref_pup_u[5] = {%s};\nconst double ref_plo_u[5] = {%s};\n", f17(pu), f17(pl)))

# Hodges-Lehmann bounds: rbounds searches a grid of step pr from the
# wilcox.test estimate
h1 <- hlsens(d, rep(0, 60), pr = 0.1, Gamma = 3, GammaInc = 0.5)$bounds
h4 <- hlsens(d, rep(0, 60), pr = 1e-4, Gamma = 3, GammaInc = 0.5)$bounds
cat("\n// hlsens(d, 0, pr = 0.1, Gamma = 3, GammaInc = 0.5)$bounds: lower, upper\n")
cat(sprintf("const double ref_hl1_lo[5] = {%s};\nconst double ref_hl1_hi[5] = {%s};\n", f17(h1[, 2]), f17(h1[, 3])))
cat("// the same with pr = 1e-4\n")
cat(sprintf("const double ref_hl4_lo[5] = {%s};\nconst double ref_hl4_hi[5] = {%s};\n", f17(h4[, 2]), f17(h4[, 3])))

# exact bounds from the Walsh averages (Rosenbaum 1993): with M = n(n+1)/2
# and T(tau) = sum of ranks of positive d - tau, the upper bound is
# (sup{tau: T > mu} + inf{tau: T < mu}) / 2 at mu = M / (1 + Gamma), the
# lower bound the same at mu = Gamma M / (1 + Gamma)
W <- sort(outer(d, d, "+")[upper.tri(diag(60), diag = TRUE)] / 2); M <- length(W)
ttau <- function(tau) { z <- d - tau; z <- z[z != 0]; sum(rank(abs(z))[z > 0]) }
ex <- function(mu) {
  a <- M - mu
  (W[ceiling(a)] + W[floor(a) + 1]) / 2
}
elo <- sapply(g, function(G) ex(G * M / (1 + G))); ehi <- sapply(g, function(G) ex(M / (1 + G)))
cat("// exact by sorting all Walsh averages in base R\n")
cat(sprintf("const double ref_ex_lo[5] = {%s};\nconst double ref_ex_hi[5] = {%s};\n", f17(elo), f17(ehi)))
# check of the Walsh identity T(tau) = #{W > tau} + #{W = tau} / 2 on a grid
taus <- seq(-3, 4, by = 1 / 64)
dev <- max(abs(sapply(taus, ttau) - sapply(taus, function(t) sum(W > t) + sum(W == t) / 2)))
cat("# max |T(tau) - (#{W > tau} + #{W = tau}/2)| over", length(taus), "grid points:", dev, "\n")
cat("# hlsens pr = 0.1 off the exact bounds by up to", max(abs(c(h1[, 2] - elo, h1[, 3] - ehi))), "\n")

set.seed(1); n <- 50000; x <- rnorm(n, 0.3); y <- rnorm(n)
cat(sprintf("\nTIME psens 50,000 pairs: %.3f s (1 run)\n",
            system.time(psens(x, y, Gamma = 2, GammaInc = 0.25))[["elapsed"]]))
cat(sprintf("TIME hlsens 50,000 pairs: %.3f s (1 run)\n",
            system.time(hlsens(x, y, Gamma = 2, GammaInc = 0.25))[["elapsed"]]))
