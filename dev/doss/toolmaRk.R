# Reference values and defect reproducer for src/lws.cpp (clab_lwss) and
# src/chm.cpp (clab_chmb). Reference implementations: stats::lowess (R
# 4.5.0) and toolmaRk 0.0.2 (CRAN).
# driver: optional path of a compiled crimlab driver (dev/doss/chmdrv.cpp)
# to compare clab_chmb with chumbley_non_random on the ameslab profiles.
source("dev/doss/pre.R")
suppressMessages(library(toolmaRk))
cat("toolmaRk", format(packageVersion("toolmaRk")), " R", format(getRversion()), "\n")

# 1. lowess cases, rebuilt in src/t_lws.cpp: case k = 0..39
ns <- c(2, 3, 4, 7, 10, 25, 60, 150); fs <- c(0.05, 0.3, 2 / 3, 1)
its <- c(0, 1, 3); ds <- c(0, 0.01, 0.2)
lcase <- function(k) {
  n <- ns[k %% 8 + 1]
  u <- minstd(2 * n, 1000 + k)
  x <- cumsum(floor(u[1:n] * 3)); i <- 0:(n - 1)
  y <- 0.5 * x + (u[n + 1:n] - 0.5) * 4 + ifelse(i %% 17 == 5, 25, 0)
  list(x = x, y = y, f = fs[(k %/% 2) %% 4 + 1], iter = its[(k %/% 3) %% 3 + 1],
       delta = ds[(k %/% 5) %% 3 + 1] * diff(range(x)))
}
out <- NULL
for (k in 0:39) { z <- lcase(k); out <- c(out, lowess(z$x, z$y, z$f, z$iter, z$delta)$y) }
cat("\n// stats::lowess(x, y, f, iter, delta) for the 40 cases, concatenated\n")
cat(sprintf("const double lw_ref[%d] = {%s};\n", length(out), f17(out)))

# 2. Chumbley score on deterministic marks, rebuilt in src/t_chm.cpp: a
# random walk w of 2,600 steps; mark 1 = w[1:1500] + noise, mark 2 =
# w[301:1800] + noise (same source, shifted 300), mark 3 = an independent
# walk (different source)
u <- minstd(3 * 2600, 4242)
w1 <- cumsum(u[1:2600] - 0.5); w3 <- cumsum(u[2601:5200] - 0.5); e <- (u[5201:7800] - 0.5) * 0.2
m1 <- w1[1:1500] + e[1:1500]; m2 <- w1[301:1800] + e[1001:2500]; m3 <- w3[1:1500] + e[1:1500]
r12 <- chumbley_non_random(matrix(m1), matrix(m2), window_opt = 300, window_val = 30, coarse = 0.25)
r13 <- chumbley_non_random(matrix(m1), matrix(m3), window_opt = 300, window_val = 30, coarse = 0.25)
r12d <- chumbley_non_random(matrix(m1), matrix(m2))
cat("\n// chumbley_non_random(matrix(m1), matrix(m2 or m3), 300, 30, 0.25): n_same, n_diff, U, p\n")
cat(sprintf("const double ch12[4] = {%s};\n", f17(unlist(r12))))
cat(sprintf("const double ch13[4] = {%s};\n", f17(unlist(r13))))
cat(sprintf("const double ch12d[4] = {%s}; // defaults 500, 50, 0.25\n", f17(unlist(r12d))))

# 3. defects (dossier toolmaRk.md section 3)
cat("\n# profile shorter than window_opt:", tryCatch({chumbley_non_random(matrix(m1[1:300]), matrix(m2[1:300])); "ran"},
    error = conditionMessage), "\n")
mn <- m1; mn[700] <- NA
cat("# NA inside a profile:", tryCatch(paste(unlist(chumbley_non_random(matrix(mn), matrix(m2), 300, 30)), collapse = " "),
    error = conditionMessage), "\n")
set.seed(1); pv <- replicate(20, chumbley_non_random(matrix(cumsum(rnorm(2000))), matrix(cumsum(rnorm(2000))))$p_value)
cat("# null p-values, 20 independent pairs: NA count", sum(is.na(pv)), "\n")

# 4. ameslab: the package's own data, through a compiled driver if given
data(ameslab)
pr <- lapply(ameslab$profile, function(d) d[[1]])
args <- commandArgs(TRUE)
pairs <- rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(1, 3), c(2, 9))
for (k in seq_len(nrow(pairs))) {
  a <- pr[[pairs[k, 1]]]; b <- pr[[pairs[k, 2]]]
  t0 <- system.time(r <- unlist(chumbley_non_random(matrix(a), matrix(b))))[["elapsed"]]
  line <- sprintf("ameslab %d-%d: toolmaRk %s (%.2f s)", pairs[k, 1], pairs[k, 2], f17(r), t0)
  if (length(args)) {
    fa <- tempfile(); fb <- tempfile(); writeLines(sprintf("%.17g", a), fa); writeLines(sprintf("%.17g", b), fb)
    cc <- system2(args[1], c(fa, fb, 500, 50, 0.25), stdout = TRUE)
    line <- paste(line, "| crimlab", cc)
    unlink(c(fa, fb))
  }
  cat(line, "\n")
}
set.seed(1); s <- cumsum(rnorm(10000))
t <- replicate(3, system.time(chumbley_non_random(matrix(s[1:9500]), matrix(s[501:10000])))[["elapsed"]])
cat(sprintf("TIME chumbley_non_random 9,500-point profiles: median %.3f s (k=3)\n", median(t)))
t <- replicate(5, system.time(chumbley_non_random(matrix(s[1:2000]), matrix(s[101:2100])))[["elapsed"]])
cat(sprintf("TIME chumbley_non_random 2,000-point profiles: median %.3f s (k=5)\n", median(t)))
