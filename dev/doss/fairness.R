# Reference values for src/fair.cpp (clab_fair).
# Field package: fairness 1.2.3 (CRAN), run on R 4.5.0.
# Prints C arrays for src/t_fair.cpp and a timing.
source("dev/doss/pre.R")
suppressMessages(library(fairness))
cat("fairness", format(packageVersion("fairness")), "\n")

# the deterministic sample src/t_fair.cpp rebuilds (i = 0, ..., 999)
i <- 0:999
y <- as.integer(((i * 13) %% 7) < 3)
g <- i %% 3
p <- (((i * 37) %% 101) + 30 * y) / 131
d <- data.frame(y = factor(y, levels = 0:1), g = c("a", "b", "c")[g + 1], p = p)

met <- list(
  acc = acc_parity, tpr = equal_odds, tnr = spec_parity, fpr = fpr_parity,
  fnr = fnr_parity, ppv = pred_rate_parity, npv = npv_parity,
  ppr = prop_parity, npos = dem_parity, mcc = mcc_parity
)
cat("\n// parity: fairness::<fn>(d, \"y\", \"g\", probs = \"p\", cutoff = 0.5)$Metric[1, ]\n")
for (k in names(met)) {
  v <- suppressWarnings(met[[k]](d, "y", "g", probs = "p", cutoff = 0.5)$Metric[1, ])
  cat(sprintf("static const double ref_%s[3] = {%s};\n", k, f17(v)))
}

# Wilson score intervals (base R stats::prop.test without continuity
# correction), the reference for clab_fair's intervals
cat("\n// prop.test(x, m, correct = FALSE, conf.level = 0.95)$conf.int\n")
xs <- c(0, 1, 7, 50, 99, 100); m <- 100
ci <- sapply(xs, function(x) suppressWarnings(prop.test(x, m, correct = FALSE)$conf.int))
cat(sprintf("static const double ref_wlo[6] = {%s};\n", f17(ci[1, ])))
cat(sprintf("static const double ref_whi[6] = {%s};\n", f17(ci[2, ])))
cat(sprintf("static const double ref_z = %s; // qnorm(0.975)\n", f17(qnorm(0.975))))

# timing
set.seed(1)
big <- data.frame(y = factor(rbinom(1e5, 1, .4), levels = 0:1),
                  g = sample(letters[1:4], 1e5, TRUE), p = runif(1e5))
t <- replicate(5, system.time(equal_odds(big, "y", "g", probs = "p"))[["elapsed"]])
cat(sprintf("\nTIME equal_odds n=1e5 4 groups: median %.3f s (k=5)\n", median(t)))
