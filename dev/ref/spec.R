# Reference values from spdep's local G*: for tests/testthat/test-spec.R,
# localG() with the cell in its own queen neighbourhood on the grid of
# counts the test builds, written to tests/testthat/ref/spec.txt; and for
# src/t_gst.cpp (clab_gstr), printed as C++ arrays, localG() on a
# deterministic count grid (spdep is the engine behind sfhotspot's
# hotspot_gistar). From the package root, with spdep in the library:
#   Rscript dev/ref/spec.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- "spdep"
source(file.path("dev", "ref", "put.R"))
source(file.path("dev", "ref", "cpp", "pre.R"))

# test-spec.R: the Gi* of gstr() against spdep's
set.seed(28)
m <- matrix(rpois(48, 3), 6, 8)
nb <- spdep::include.self(spdep::cell2nb(6, 8, type = "queen"))
lg <- spdep::localG(as.vector(t(m)), spdep::nb2listw(nb, style = "B"))

# src/t_gst.cpp: 12 rows x 15 columns, cell (r, c) at centre (c, r),
# counts column-major (r fastest): floor(12 u) plus 20 in a 3 x 3 block
# centred at (r, c) = (4, 6)
cat("spdep", format(packageVersion("spdep")), "\n")
nr <- 12
nc <- 15
u <- minstd(nr * nc, 99)
cnt <- floor(12 * u)
rr <- rep(0:(nr - 1), nc)
cc <- rep(0:(nc - 1), each = nr)
cnt <- cnt + ifelse(abs(rr - 4) <= 1 & abs(cc - 6) <= 1, 20, 0)
xy <- cbind(cc, rr)
for (d in c(1.5, 2.5)) {
  g <- as.numeric(spdep::localG(cnt, spdep::nb2listw(
    spdep::include.self(spdep::dnearneigh(xy, 0, d)), style = "B")))
  cat(sprintf(paste0("// localG(cnt, nb2listw(include.self(dnearneigh(xy, ",
                     "0, %g)), style = \"B\"))\n"), d))
  cat(sprintf("const double gref%d[%d] = {%s};\n", round(10 * d), length(g),
              f17(g)))
}

put(list(inp = sum(m), lg = as.vector(lg)), "spec", need)
