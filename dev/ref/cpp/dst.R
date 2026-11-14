# Reference values for src/t_nnd.cpp and src/t_rpk.cpp (nearest-neighbour
# index, k-order distances, Ripley's K). Run on zeus:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/ref/cpp/dst.R
suppressMessages(library(spatstat.explore))
g <- function(i, a, m) ((i * a) %% m) / m
pr <- function(name, v) cat(name, " = {", paste(sprintf("%.17g", v),
  collapse = ", "), "};\n", sep = "")
cat("// R", as.character(getRversion()), "spatstat.explore",
  as.character(packageVersion("spatstat.explore")), "spatstat.geom",
  as.character(packageVersion("spatstat.geom")), "\n")
# pattern P: 150 points in the rectangle (0, 2) x (0, 1) and 5 copies of
# its first 5 points
i <- 0:149
# a quadratic recurrence (not a lattice); exact in doubles below 2^53
x <- ((i * i * 7919 + i * 31) %% 1009) / 1009 * 2
y <- ((i * i * 104729 + i * 7) %% 1013) / 1013
W <- owin(c(0, 2), c(0, 1))
XD <- ppp(c(x, x[1:5]), c(y, y[1:5]), window = W, check = FALSE)
cat("// colMeans(nndist(XD, k = 1:4))\n")
pr("nnmean", colMeans(nndist(XD, k = 1:4)))
cat("// nndist(XD, k = 1:4)[1:10, ] by row\n")
pr("nnrows", t(nndist(XD, k = 1:4)[1:10, ]))
X <- ppp(x, y, window = W)
cat("// clarkevans(X, correction = c(\"none\", \"Donnelly\"))\n")
pr("ce", clarkevans(X, correction = c("none", "Donnelly")))
r <- seq(0, 0.25, by = 0.005)
k <- as.data.frame(Kest(X, r = r,
  correction = c("none", "border", "bord.modif", "isotropic")))
cat("// Kest(X, r = seq(0, 0.25, by = 0.005), correction = ...)\n")
pr("kun", k$un)
pr("kbord", k$border)
# the modified border estimator from spatstat's dedicated path: asked for
# together with "isotropic" or "none", Kest divides by the window's area
# where the definition (and the dedicated path) divides by the eroded area
kbm <- Kest(X, r = r, correction = "bord.modif")$bord.modif
pr("kbm", kbm)
kbm2 <- k$bord.modif
cat("// Kest(..., correction = c(\"none\", ..., \"isotropic\"))$bord.modif",
  "times |W| / |W eroded by r| equals the dedicated value to",
  sprintf("%.1e", max(abs(kbm2 * 2 / ((2 - 2 * r) * (1 - 2 * r)) - kbm)[r < 0.5])),
  "\n")
pr("kiso", k$iso)
