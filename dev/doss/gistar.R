# Reference values for src/gst.cpp (clab_gstr): spdep's local G* on a
# deterministic count grid (spdep is the engine behind sfhotspot's
# hotspot_gistar), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/gistar.R
suppressMessages(library(spdep))
cat("spdep", format(packageVersion("spdep")), "\n")
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
minstd <- function(n, s) { u <- numeric(n); for (i in 1:n) { s <- (48271 * s) %% 2147483647; u[i] <- s / 2147483647 }; u }
# 12 rows x 15 columns, cell (r, c) at centre (c, r), counts column-major
# (r fastest): floor(12 u) plus 20 in a 3 x 3 block centred at (r, c) = (4, 6)
nr <- 12; nc <- 15
u <- minstd(nr * nc, 99)
cnt <- floor(12 * u)
rr <- rep(0:(nr - 1), nc); cc <- rep(0:(nc - 1), each = nr)
cnt <- cnt + ifelse(abs(rr - 4) <= 1 & abs(cc - 6) <= 1, 20, 0)
xy <- cbind(cc, rr)
for (d in c(1.5, 2.5)) {
  nb <- include.self(dnearneigh(xy, 0, d))
  g <- as.numeric(localG(cnt, nb2listw(nb, style = "B")))
  cat(sprintf("// localG(cnt, nb2listw(include.self(dnearneigh(xy, 0, %g)), style = \"B\"))\n", d))
  cat(sprintf("const double gref%d[%d] = {%s};\n", round(10 * d), length(g), f17(g)))
}
