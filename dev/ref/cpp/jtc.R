# Reference values for src/t_jtc.cpp (decay functions, trip-length kernel
# densities). Base R only. Run on zeus:
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null nice -n 10 \
#   /usr/bin/Rscript dev/ref/cpp/jtc.R
pr <- function(name, v) cat(name, " = {", paste(sprintf("%.17g", v),
  collapse = ", "), "};\n", sep = "")
cat("// R", as.character(getRversion()), "\n")
d <- c(0, 0.1, 0.5, 1, 2.5, 4.2, 7, 12.5, 30)
cat("// 2.5 * dnorm(d, 4.2, 4.6)\n")
pr("rnorm", 2.5 * dnorm(d, 4.2, 4.6))
cat("// 8.6 * dlnorm(d, 0.7, 0.9)\n")
pr("rlnorm", 8.6 * dlnorm(d, 0.7, 0.9))
cat("// 1.89 * dexp(d, 0.06) / 0.06\n")
pr("rnexp", 1.89 * dexp(d, 0.06) / 0.06)
# trip lengths and kernel densities written out from their definitions
i <- 0:59
tr <- ((i * i * 31 + 7 * i) %% 97) / 10
t <- seq(0, 10, by = 0.5)
h <- 0.8
cat("// trips tr = ((i^2 31 + 7 i) mod 97) / 10, i = 0..59; t = 0, 0.5, ..., 10; h = 0.8\n")
cat("// normal kernel: sapply(t, function(s) sum(dnorm((s - tr) / h)) / (60 h))\n")
pr("kdn", sapply(t, function(s) sum(dnorm((s - tr) / h)) / (60 * h)))
cat("// reflected: adds dnorm((s + tr) / h)\n")
pr("kdnr", sapply(t, function(s) sum(dnorm((s - tr) / h) + dnorm((s + tr) / h)) / (60 * h)))
q <- function(u) ifelse(abs(u) < 1, 15 / 16 * (1 - u^2)^2, 0)
cat("// quartic kernel 15/16 (1 - u^2)^2\n")
pr("kdq", sapply(t, function(s) sum(q((s - tr) / h)) / (60 * h)))
