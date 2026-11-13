# Test error of rfor() settings and grf, like for like, over data-generating
# functions and seeds:  Rscript dev/rftune.R  (writes dev/rftune.csv)
suppressPackageStartupMessages(library(crimlab))
stopifnot(requireNamespace("grf"))
ntree <- 200
nseed <- as.integer(Sys.getenv("RFT_SEEDS", "10"))
gen <- list(
  sparse = function(x) drop(x[, 1:4] %*% c(1, -0.5, 0.5, 0.25)),
  friedman = function(x) 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5],
  step = function(x) 2 * sin(x[, 1]) + (x[, 2] > 0) - 0.5 * x[, 3]^2 + x[, 4],
  prob = function(x) plogis(x[, 1] - x[, 2] + 0.5 * x[, 3] * x[, 4])
)
cfg <- list(
  # name, rfor arguments, grf arguments (NULL: no grf twin)
  list("rf-default", list(repl = TRUE, frac = 1, minnode = 5, hon = FALSE), NULL, "sqrt"),
  list("adaptive-half", list(repl = FALSE, frac = 0.5, minnode = 5, hon = FALSE, imb = 0.05),
       list(honesty = FALSE), "grf"),
  list("honest-half", list(repl = FALSE, frac = 0.5, minnode = 5, hon = TRUE, imb = 0.05),
       list(honesty = TRUE), "grf"),
  list("adaptive-half-sqrt", list(repl = FALSE, frac = 0.5, minnode = 5, hon = FALSE, imb = 0.05),
       list(honesty = FALSE), "sqrt"),
  list("rf-grfmtry", list(repl = TRUE, frac = 1, minnode = 5, hon = FALSE), NULL, "grf")
)
out <- list()
for (p in c(20, 200)) for (g in names(gen)) for (s in seq_len(nseed)) {
  set.seed(1000 * s + p)
  n <- 2000
  x <- matrix(if (g == "friedman") runif(2 * n * p) else rnorm(2 * n * p), 2 * n)
  f <- gen[[g]](x)
  y <- if (g == "prob") rbinom(2 * n, 1, f) else f + rnorm(2 * n)
  tr <- 1:n; te <- n + 1:n
  for (c in cfg) {
    mt <- if (c[[4]] == "sqrt") floor(sqrt(p)) else min(ceiling(sqrt(p) + 20), p)
    a <- do.call(rfor, c(list(x[tr, ], y[tr], x[te, ], ntree = ntree, mtry = mt, seed = s), c[[2]]))
    out[[length(out) + 1]] <- data.frame(p = p, gen = g, seed = s, cfg = c[[1]], pkg = "crimlab",
      mse = mean((a$pred - f[te])^2))
    if (!is.null(c[[3]])) {
      r <- do.call(grf::regression_forest, c(list(x[tr, ], y[tr], num.trees = ntree, mtry = mt,
        min.node.size = 5, sample.fraction = 0.5, alpha = 0.05, ci.group.size = 1,
        num.threads = 1, seed = s), c[[3]]))
      out[[length(out) + 1]] <- data.frame(p = p, gen = g, seed = s, cfg = c[[1]], pkg = "grf",
        mse = mean((predict(r, x[te, ])$predictions - f[te])^2))
    }
  }
}
d <- do.call(rbind, out)
utils::write.csv(d, "dev/rftune.csv", row.names = FALSE)
a <- aggregate(mse ~ p + gen + cfg + pkg, d, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
print(a[order(a$p, a$gen, a$cfg, a$pkg), ], digits = 3)
