# Reference values for src/tfr.cpp (clab_tfrp). Field package: tfer 1.3
# (CRAN), R 4.5.0.
source("dev/doss/pre.R")
library(tfer)
cat("tfer", format(packageVersion("tfer")), "\n")
# reference distribution for src/t_tfr.cpp: the default model, 10^6 draws
# (tfer's own Monte Carlo; standard error sqrt(p (1 - p) / 1e6))
set.seed(2026)
big <- transfer(N = 1e6)
p <- tprob(big)
cat("\n// tfer 1.3: set.seed(2026); tprob(transfer(N = 1e6))[1:81] (counts 0..80)\n")
cat(sprintf("const double tf_ref[81] = {%s};\n", paste(sprintf("%.17g", p[1:81]), collapse = ", ")))
cat("// mean of the 1e6 simulated counts:", mean(getValues(big)), "\n")
cat("# max simulated count:", length(p) - 1, "\n")
t <- replicate(5, system.time(transfer(N = 1e5))[["elapsed"]])
cat(sprintf("TIME transfer N=1e5: median %.3f s (k=5)\n", median(t)))
