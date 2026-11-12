# Defect reproducer for the transfer model (no crimlab kernel yet; see
# dev/wip-h-doss.md). Field package: tfer 1.3 (CRAN), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/tfer.R
library(tfer)
cat("tfer", format(packageVersion("tfer")), "\n")
set.seed(1); a <- transfer(N = 2000); m <- length(tprob(a)) - 1
cat("max simulated count m =", m, "\n")
cat("tprob(a, c(2, m + 50)):\n"); print(tprob(a, c(2, m + 50)))
cat("tprob(a, c(5, 1)) (input order not kept):\n"); print(tprob(a, c(5, 1)))
cat("Q = -1:", tryCatch(format(transfer(N = 10, Q = -1)$results), error = conditionMessage,
    warning = conditionMessage), "\n")
# the per-step loss rates are drawn once per case and shared by all its
# fragments, so the recovered count is a Poisson mixture, not a Poisson
# thinning: overdispersion against a Poisson with the same mean
v <- getValues(transfer(N = 20000))
cat("mean", mean(v), "variance", var(v), "share of zeros", mean(v == 0), "\n")
