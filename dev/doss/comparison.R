# Reference values and defect reproducer for src/tlv.cpp (clab_tlvc,
# clab_tlvlr). Field package: comparison 1.0.8 (CRAN), R 4.5.0.
source("dev/doss/pre.R")
suppressMessages(library(comparison))
cat("comparison", format(packageVersion("comparison")), "\n")

# 1. the documented glass example (K, Ca, Fe)
data(glass)
Z <- two.level.components(glass, c(7, 8, 9), 1)
ctl <- two.level.comparison.items(glass[1:6, ], c(7, 8, 9))
r1 <- two.level.comparison.items(glass[7:12, ], c(7, 8, 9))
r2 <- two.level.comparison.items(glass[19:24, ], c(7, 8, 9))
cat("\n// comparison 1.0.8: Z <- two.level.components(glass, c(7, 8, 9), 1);\n")
cat("// items from glass rows 1:6 (control), 7:12 (r1), 19:24 (r2), n = 6\n")
carr("gU", Z$v.within); carr("gC", Z$v.between); carr("gmu", Z$overall.means)
carr("gctl", ctl$item.means); carr("gr1", r1$item.means); carr("gr2", r2$item.means)
cat("// two.level.normal.LR(ctl, r1, Z), (ctl, r2, Z)\n")
carr("glr", c(two.level.normal.LR(ctl, r1, Z), two.level.normal.LR(ctl, r2, Z)))

# 2. defect: natural-scale LR (dossier comparison.md section 3)
sh <- function(it, k) { it$item.means <- it$item.means + k * sqrt(diag(Z$v.between)); it }
cat("\n# defect: both items 40 between-SD from the mean (same source): ",
    two.level.normal.LR(sh(ctl, 40), sh(ctl, 40), Z), "\n")
cat("# defect: items 40 between-SD apart: ",
    two.level.normal.LR(sh(ctl, 20), sh(ctl, -20), Z), "\n")
# an independent log-space value for the same inputs: the joint density of
# the two item means under each hypothesis (2p-variate normal)
ljoint <- function(a, b, U, C, mu, na, nb, same) {
  p <- length(a)
  S <- rbind(cbind(U / na + C, if (same) C else 0 * C),
             cbind(if (same) C else 0 * C, U / nb + C))
  r <- c(a - mu, b - mu)
  -0.5 * (determinant(S)$modulus + sum(r * solve(S, r)) + 2 * p * log(2 * pi))
}
llr <- function(a, b) ljoint(a, b, Z$v.within, Z$v.between, Z$overall.means, 6, 6, TRUE) -
  ljoint(a, b, Z$v.within, Z$v.between, Z$overall.means, 6, 6, FALSE)
s40 <- sh(ctl, 40)$item.means; p20 <- sh(ctl, 20)$item.means; m20 <- sh(ctl, -20)$item.means
cat("// natural log LR by the 2p-variate joint density (base R solve/determinant)\n")
carr("gfar", c(llr(s40, s40), llr(p20, m20)))
carr("gs40", s40); carr("gp20", p20); carr("gm20", m20)
cat("// check of that route on the finite example: log(glr) =", f17(log(two.level.normal.LR(ctl, r1, Z))),
    " joint:", f17(llr(ctl$item.means, r1$item.means)), "\n")

# 3. components on deterministic data (MINSTD generator, rebuilt in C)
mk <- function(m, nrep) {
  u <- minstd(sum(nrep) * 3 + m * 3, 12345)
  eff <- matrix(u[1:(m * 3)], m, 3) * c(1, 2, 0.5)
  rows <- NULL; k <- m * 3
  for (i in 1:m) for (j in 1:nrep[i]) {
    rows <- rbind(rows, c(i, eff[i, ] + 0.1 * (u[k + 1:3] - 0.5))); k <- k + 3
  }
  d <- as.data.frame(rows); names(d) <- c("item", "a", "b", "c"); d$item <- factor(d$item); d
}
bal <- mk(30, rep(4, 30))
Zb <- two.level.components(bal, 2:4, 1)
cat("\n// balanced: m = 30 items x 4, two.level.components(d, 2:4, 1)\n")
carr("bU", Zb$v.within); carr("bC", Zb$v.between); carr("bmu", Zb$overall.means)
nub <- rep(c(2, 3, 5, 8), length.out = 30)
ub <- mk(30, nub)
Zu <- two.level.components(ub, 2:4, 1)
cat("// unbalanced: n_i = 2, 3, 5, 8 repeating; what comparison reports\n")
carr("uU", Zu$v.within); carr("uC", Zu$v.between)
# the standard (Searle et al. 1992) estimates for the unbalanced data
X <- as.matrix(ub[, 2:4]); g <- ub$item; N <- nrow(X); m <- nlevels(g)
mi <- rowsum(X, g) / as.vector(table(g)); ni <- as.vector(table(g))
W <- crossprod(X - mi[as.integer(g), ]); U <- W / (N - m)
gm <- colMeans(X); SB <- crossprod((mi - rep(gm, each = m)) * sqrt(ni))
C <- (SB - (m - 1) * U) / (N - sum(ni^2) / N)
cat("// ANOVA estimates by base R (crossprod, rowsum)\n")
carr("sU", U); carr("sC", C)

# 4. bias of comparison's within estimate on unbalanced data, by simulation:
# true U = diag(1, 1, 1), C = diag(4, 4, 4); 400 data sets
set.seed(7); acc <- matrix(0, 3, 3); acc2 <- acc
for (r in 1:400) {
  th <- matrix(rnorm(30 * 3, 0, 2), 30, 3)
  X <- th[rep(1:30, nub), ] + matrix(rnorm(sum(nub) * 3), ncol = 3)
  d <- data.frame(item = factor(rep(1:30, nub)), X)
  acc <- acc + suppressWarnings(two.level.components(d, 2:4, 1))$v.within / 400
  mi <- rowsum(X, d$item) / nub
  acc2 <- acc2 + crossprod(X - mi[rep(1:30, nub), ]) / (sum(nub) - 30) / 400
}
cat("\n# unbalanced bias, mean of 400 within estimates (truth 1): comparison",
    round(diag(acc), 3), " pooled within:", round(diag(acc2), 3), "\n")
cat("# expected factor (m/N) sum n_i (n_i - 1) / (N - m) =",
    (30 / sum(nub)) * sum(nub * (nub - 1)) / (sum(nub) - 30), "\n")

t <- replicate(5, system.time(two.level.components(glass, c(7, 8, 9), 1))[["elapsed"]])
cat(sprintf("\nTIME two.level.components glass 200 items x 12 x 3: median %.3f s (k=5)\n", median(t)))
t <- replicate(5, system.time(for (k in 1:100) two.level.normal.LR(ctl, r1, Z))[["elapsed"]])
cat(sprintf("TIME two.level.normal.LR: median %.6f s per call (k=5 x 100)\n", median(t) / 100))
