# Shared by the reference generators dev/ref/cpp/*.R and dev/ref/spec.R; run each from the package root on zeus:
# Run from the package root: Rscript dev/ref/cpp/<name>.R
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
carr <- function(nm, x) cat(sprintf("const double %s[%d] = {%s};\n", nm, length(x), f17(x)))
minstd <- function(n, s) { u <- numeric(n); for (i in 1:n) { s <- (48271 * s) %% 2147483647; u[i] <- s / 2147483647 }; u }
