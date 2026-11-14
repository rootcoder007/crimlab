# Reference values for src/t_cgr.cpp and src/t_geo.cpp (centrography,
# convex hull, Weber point, great-circle distances). Run on zeus:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/ref/cpp/cgr.R
# The point sets are built from integer recurrences that the C++ tests
# reproduce operation for operation, so both sides see the same doubles.
g <- function(i, a, m) ((i * a) %% m) / m
pr <- function(name, v) cat(name, " = {", paste(sprintf("%.17g", v),
  collapse = ", "), "};\n", sep = "")
cat("// R", as.character(getRversion()), "\n")

# 1. convex hull: base R grDevices::chull (Eddy 1977), 1-based indices in
# clockwise order
i <- 0:199
x <- g(i, 7919, 1009) * 1000
y <- g(i, 104729, 1013) * 1000
h <- chull(x, y) - 1
cat("// grDevices::chull(x, y) - 1, set A (n = 200)\n")
pr("hullA", h)
# the triple (187, 5, 109): exactly collinear in the rationals, not in
# their doubles (sign of the orientation on each)
library(gmp)
k3 <- c(187, 5, 109)
or3 <- function(u, v) sign((u[2] - u[1]) * (v[3] - v[1]) - (v[2] - v[1]) * (u[3] - u[1]))
cat("// orientation of (187, 5, 109): rationals",
  as.character(or3(as.bigq((k3 * 7919) %% 1009, 1009) * 1000,
    as.bigq((k3 * 104729) %% 1013, 1013) * 1000)),
  "doubles", as.character(or3(as.bigq(x[k3 + 1]), as.bigq(y[k3 + 1]))), "\n")
# lattice 10 x 10 with every point duplicated: corners only
j <- 0:199
x2 <- (j %% 100) %% 10
y2 <- (j %% 100) %/% 10
cat("// grDevices::chull, set B (10 x 10 lattice, each point twice)\n")
pr("hullBxy", unlist(lapply(chull(x2, y2), function(k) c(x2[k], y2[k]))))

# 2. Weber point: pcaPP::l1median (Croux, Filzmoser and Fritz), default
# method, on set A with weights replicated as frequencies
library(pcaPP)
w <- 1 + (i %% 3)
X <- cbind(rep(x, w), rep(y, w))
cat("// pcaPP", as.character(packageVersion("pcaPP")),
  "l1median_NLM(cbind(rep(x, w), rep(y, w)), maxit = 1000, tol = 1e-14)$par\n")
m <- l1median_NLM(X, maxit = 1000, tol = 1e-14)$par; pr("l1A", m)
cat("// its objective sum w_i |p_i - c|\n")
pr("l1Aobj", sum(sqrt((X[, 1] - m[1])^2 + (X[, 2] - m[2])^2)))

# 3. great-circle distances: geosphere::distVincentySphere (sphere case of
# Vincenty's formula) with r = 6371008.8
library(geosphere)
k <- 0:49
p1 <- cbind(g(k, 7919, 1009) * 360 - 180, g(k, 104729, 1013) * 180 - 90)
p2 <- cbind(g(k, 31, 97) * 360 - 180, g(k, 61, 101) * 170 - 85)
cat("// geosphere", as.character(packageVersion("geosphere")),
  "distVincentySphere(p1, p2, r = 6371008.8)\n")
pr("gcd", distVincentySphere(p1, p2, r = 6371008.8))
cat("// geosphere distHaversine for the same pairs\n")
pr("gch", distHaversine(p1, p2, r = 6371008.8))
