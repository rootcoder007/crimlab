# Reference values of tests/testthat/test-sp.R: geosphere's great-circle
# distances, sf's (PROJ's) projection, spatstat's nearest neighbours,
# Clark-Evans ratios and K, and sf point objects (one read back from a
# GeoPackage), on the inputs the tests build. From the package root, with
# crimlab, geosphere, sf, spatstat.geom and spatstat.explore in the
# library:
#   Rscript dev/ref/sp.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("geosphere", "sf", "spatstat.geom", "spatstat.explore")
source(file.path("dev", "ref", "put.R"))
library(crimlab)

# "great-circle distances agree with geosphere and closed forms"
r <- 6371008.8
set.seed(3)
a <- cbind(lon = runif(100, -180, 180), lat = runif(100, -90, 90))
b <- cbind(lon = runif(100, -180, 180), lat = runif(100, -90, 90))
gc <- list(
  inp = sum(a, b),
  g = geosphere::distVincentySphere(a, b, r = r),
  g1 = geosphere::distVincentySphere(a, b[1, ], r = r)
)

# "the projection inverts and agrees with PROJ through sf"
set.seed(4)
lon <- runif(200, -100, -60)
lat <- runif(200, 20, 60)
a <- attr(eqproj(lon, lat), "centre")
crs <- sprintf("+proj=laea +lat_0=%.17g +lon_0=%.17g +R=6371008.8",
               a[2], a[1])
s <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(
  data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326), crs))
proj <- list(inp = sum(lon, lat), x = unname(s[, 1]), y = unname(s[, 2]))

# "nearest neighbours and K agree with spatstat"
set.seed(8)
x <- runif(150)
y <- 2 * runif(150)
rr <- seq(0, 0.3, by = 0.01)
X <- spatstat.geom::ppp(x, y, c(0, 1), c(0, 2))
ce <- spatstat.explore::clarkevans(X, correction = c("none", "Donnelly"))
ke <- function(co) {
  as.data.frame(spatstat.explore::Kest(X, r = rr, correction = co))
}
nn <- list(
  inp = sum(x, y),
  dist = unname(spatstat.geom::nndist(X, k = 1:3)),
  naive = ce[["naive"]],
  donnelly = ce[["Donnelly"]],
  iso = ke("isotropic")$iso,
  border = ke("border")$border,
  none = ke("none")$un,
  modif = ke("bord.modif")$bord.modif
)

# "sf points go in, with their CRS checked against lonlat"
set.seed(10)
d <- data.frame(lon = runif(20, -80, -79), lat = runif(20, 43, 44))
s <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
f <- tempfile(fileext = ".gpkg")
sf::st_write(s, f, quiet = TRUE)
g <- sf::st_read(f, quiet = TRUE)
unlink(f)
p <- sf::st_transform(s, 32617)
e <- sf::st_as_sf(eqproj(d), coords = c("x", "y"))
obj <- list(
  inp = sum(d),
  s = s,
  g = g,
  p = p,
  pg = sf::st_geometry(p),
  xy = sf::st_coordinates(p),
  ex = sf::st_coordinates(e)[, "X"],
  u = sf::st_sfc(sf::st_point(c(1, 2, 3)), sf::st_point(c(1, 2, 4))),
  l = sf::st_as_sf(sf::st_sfc(sf::st_linestring(matrix(1:4, 2))))
)

put(list(gc = gc, proj = proj, nn = nn, obj = obj), "sp", need)
