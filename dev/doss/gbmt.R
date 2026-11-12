# Reference values for src/mvt.cpp (clab_mvtll,
# clab_mvtem). Field package: gbmt 0.1.4 (CRAN), R 4.5.0.
source("dev/doss/pre.R")
suppressMessages(library(gbmt))
cat("gbmt", format(packageVersion("gbmt")), "\n")

# 60 units x 6 times x 2 indicators, two groups; values rounded to 1/8 so
# that they print exactly
set.seed(5)
n <- 60; nt <- 6; g <- rep(1:2, c(25, 35))
mk <- function(i, t) {
  if (g[i] == 1) c(1 + 0.5 * t, 3 - 0.2 * t) else c(4 - 0.3 * t, 1 + 0.1 * t)
}
d <- do.call(rbind, lapply(1:n, function(i) t(sapply(1:nt, function(t)
  c(i, t, round((mk(i, t) + rnorm(2, 0, c(0.6, 0.4))) * 8) / 8)))))
d <- as.data.frame(d); names(d) <- c("id", "time", "y1", "y2")
cat("\n// d: rows (id, time, y1, y2), id-major\n")
carr("D", t(as.matrix(d)))
fit <- gbmt(c("y1", "y2"), "id", "time", ng = 2, d = 1, data = d, scaling = 0,
            pruning = FALSE, quiet = TRUE)
cat("// gbmt(c(\"y1\", \"y2\"), \"id\", \"time\", ng = 2, d = 1, data = d, scaling = 0, pruning = FALSE)\n")
for (k in 1:2) {
  carr(sprintf("gbeta%d", k), fit$beta[[k]])      # 2 x 2: rows (intercept, t), cols (y1, y2)
  carr(sprintf("gS%d", k), fit$Sigma[[k]])
}
carr("gprob", colMeans(fit$posterior))
carr("gpost", fit$posterior[1:5, ])
cat(sprintf("const double gll = %s;\n", f17(fit$logLik)))
cat("# components:", paste(names(fit), collapse = " "), "\n")

# timing
set.seed(1)
big <- function(n, ng) {
  gg <- sample(ng, n, TRUE)
  dd <- do.call(rbind, lapply(1:n, function(i) t(sapply(1:8, function(t)
    c(i, t, c(gg[i], -gg[i]) * 0.5 * t + rnorm(2))))))
  dd <- as.data.frame(dd); names(dd) <- c("id", "time", "y1", "y2"); dd
}
b1 <- big(200, 2)
cat(sprintf("\nTIME gbmt 200 units ng=2 d=1: %.3f s (1 run)\n",
    system.time(gbmt(c("y1", "y2"), "id", "time", ng = 2, d = 1, data = b1, quiet = TRUE))[["elapsed"]]))
