# Reference values and defect reproducer for src/fair.cpp (clab_fair).
# Field package: fairness 1.2.3 (CRAN), run on R 4.5.0.
# Prints: the defect, C arrays for src/t_fair.cpp, and a timing.
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

# defect (dossier fairness.md section 3): every score above the cutoff
set.seed(1); n <- 1000
dp <- data.frame(y = factor(rbinom(n, 1, .4), levels = 0:1),
                 g = sample(c("a", "b"), n, TRUE), p = runif(n, 0.6, 1))
cat("\n# defect: all scores > cutoff 0.5\n")
cat("equal_odds sensitivity:", equal_odds(dp, "y", "g", probs = "p")$Metric[1, ], "\n")
cat("truth (hand):", tapply(dp$p > .5 & dp$y == 1, dp$g, sum) /
  tapply(dp$y == 1, dp$g, sum), "\n")
cat("dem_parity positives:", dem_parity(dp, "y", "g", probs = "p")$Metric[1, ], "\n")
cat("truth (hand):", tapply(dp$p > .5, dp$g, sum), "\n")
# the same defect on the deterministic sample, for the C test
cat("deterministic sample with cutoff -1: equal_odds",
    equal_odds(d, "y", "g", probs = "p", cutoff = -1)$Metric[1, ],
    " dem_parity", dem_parity(d, "y", "g", probs = "p", cutoff = -1)$Metric[1, ], "\n")

# timing: the dossier's standard case
set.seed(1)
big <- data.frame(y = factor(rbinom(1e5, 1, .4), levels = 0:1),
                  g = sample(letters[1:4], 1e5, TRUE), p = runif(1e5))
t <- replicate(5, system.time(equal_odds(big, "y", "g", probs = "p"))[["elapsed"]])
cat(sprintf("\nTIME equal_odds n=1e5 4 groups: median %.3f s (k=5)\n", median(t)))
