# Defect probe for Gi* calibration (no crimlab kernel yet; see
# dev/wip-h-doss.md). Field package: sfhotspot 1.1.1 (CRAN), R 4.5.0.
# Under complete spatial randomness on [0, 1000]^2 with 50 m cells, the
# share of cells with unadjusted p < 0.05, overall and for core cells whose
# 3 x 3 neighbourhood lies inside the study square.
source("dev/doss/pre.R")
suppressMessages({library(sfhotspot); library(sf)})
cat("sfhotspot", format(packageVersion("sfhotspot")), "\n")
pts <- function(n, seed) { set.seed(seed); st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)), coords = c("x", "y"), crs = 27700) }
res <- NULL
for (s in 1:20) {
  r <- hotspot_gistar(pts(4000, 200 + s), cell_size = 50, kde = FALSE, p_adjust_method = "none", quiet = TRUE)
  cc <- st_coordinates(suppressWarnings(st_centroid(r)))
  core <- cc[, 1] - 75 >= 0 & cc[, 1] + 75 <= 1000 & cc[, 2] - 75 >= 0 & cc[, 2] + 75 <= 1000
  res <- rbind(res, c(all = mean(r$pvalue < 0.05), core = mean(r$pvalue[core] < 0.05),
                      ncells = nrow(r), ncore = sum(core)))
}
cat("grid bounding box of the last run:", st_bbox(r), "\n")
cat("mean share flagged at 0.05, 20 runs: all cells", mean(res[, 1]), " core cells", mean(res[, 2]),
    " (cells", res[1, 3], ", core", res[1, 4], ")\n")
