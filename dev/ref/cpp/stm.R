# Reference values for src/t_stm.cpp (Mantel correlation). Run on zeus:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/ref/cpp/stm.R
suppressMessages(library(vegan))
pr <- function(name, v) cat(name, " = {", paste(sprintf("%.17g", v),
  collapse = ", "), "};\n", sep = "")
cat("// R", as.character(getRversion()), "vegan", as.character(packageVersion("vegan")), "\n")
i <- 0:59
x <- ((i * i * 7919 + i * 31) %% 1009) / 1009 * 2
y <- ((i * i * 104729 + i * 7) %% 1013) / 1013
t <- ((i * 17) %% 61) + x  # times related to x: a space-time link
ds <- dist(cbind(x, y)); dt <- dist(t)
cat("// mantel(dist(cbind(x, y)), dist(t), permutations = 0)$statistic\n")
pr("mr", mantel(ds, dt, permutations = 0)$statistic)
cat("// the same with 1 / (d + 0.1) and 1 / (dt + 2)\n")
pr("mrr", mantel(1 / (ds + 0.1), 1 / (dt + 2), permutations = 0)$statistic)
