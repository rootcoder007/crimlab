# Reference values for src/t_hsp.cpp (K-means). Base R only. Run on zeus:
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null nice -n 10 \
#   /usr/bin/Rscript dev/doss/hsp.R
pr <- function(name, v) cat(name, " = {", paste(sprintf("%.17g", v),
  collapse = ", "), "};\n", sep = "")
cat("// R", as.character(getRversion()), "\n")
i <- 0:199
x <- ((i * i * 7919 + i * 31) %% 1009) / 1009 * 2
y <- ((i * i * 104729 + i * 7) %% 1013) / 1013
X <- cbind(x, y)
for (k in c(3, 6)) {
  km <- kmeans(X, centers = X[1:k, ], algorithm = "Lloyd", iter.max = 1000)
  cat(sprintf("// kmeans(X, centers = X[1:%d, ], algorithm = \"Lloyd\", iter.max = 1000)\n", k))
  pr(sprintf("km%dcl", k), km$cluster - 1)
  pr(sprintf("km%dcen", k), c(km$centers[, 1], km$centers[, 2]))
  pr(sprintf("km%dwss", k), c(km$tot.withinss, km$withinss))
}
