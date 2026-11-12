# Shared by the reference generators dev/doss/*.R; run each from the package root on zeus:
# R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/<pkg>.R
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
carr <- function(nm, x) cat(sprintf("const double %s[%d] = {%s};\n", nm, length(x), f17(x)))
minstd <- function(n, s) { u <- numeric(n); for (i in 1:n) { s <- (48271 * s) %% 2147483647; u[i] <- s / 2147483647 }; u }
