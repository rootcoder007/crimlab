# Reference values and defect reproducer for src/spt.cpp (clab_sptst).
# Field package: sppt 1.0.0 (CRAN), R 4.5.0.
source("dev/doss/pre.R")
suppressMessages({library(sppt); library(sf)})
cat("sppt", format(packageVersion("sppt")), "\n")
grid <- st_make_grid(st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))))), n = c(10, 10))
# deterministic counts, rebuilt in src/t_spt.cpp: a = floor(10 u), b =
# floor(10 u') for units 0..99 (MINSTD from 31), b + 30 at unit 17, both
# 0 at units 40 and 41
u <- minstd(200, 31); i <- 0:99
a <- floor(10 * u[1:100]); b <- floor(10 * u[101:200])
b[18] <- b[18] + 30; a[41:42] <- 0; b[41:42] <- 0
d <- st_sf(group = 1:100, a = a, b = b, geometry = grid)
r <- suppressWarnings(sppt(d, count_col = c("a", "b"), B = 4999, seed = 1, arrangement = FALSE,
                           create_maps = FALSE, verbose = FALSE))
s <- attr(r, "sppt_summary")
cat("\n// sppt(d, count_col = c(\"a\", \"b\"), B = 4999, seed = 1): difference, pvalue (NA -> NaN)\n")
cat(sprintf("const double sp_diff[100] = {%s};\n", f17(r$difference)))
cat(sprintf("const double sp_p[100] = {%s};\n", gsub("NA", "NAN", f17(r$pvalue))))
cat("# summary names:", names(s), "\n")
cat("# S index:", s$s_index_pv, " n_sig:", s$n_sig_pv, " eligible:", s$n_eligible,
    " sameness p:", s$sameness$p_value, s$sameness$method, "\n")
# the bootstrap p-value is unstable beyond Monte Carlo error at sparse
# units: replicates fall on the lattice of share steps, and whether one
# ties the observed difference after centring on the replicate mean is
# decided by the noise of that mean. Range of p over 8 seeds, B = 4999:
pp <- sapply(1:8, function(sd) suppressWarnings(sppt(d, count_col = c("a", "b"), B = 4999, seed = sd,
  arrangement = FALSE, create_maps = FALSE, verbose = FALSE))$pvalue)
cat("// sppt(d, ..., B = 4999, seed = 1..8): per-unit minimum and maximum p\n")
cat(sprintf("const double sp_pmin[100] = {%s};\n", gsub("NA", "NAN", f17(apply(pp, 1, min)))))
cat(sprintf("const double sp_pmax[100] = {%s};\n", gsub("NA", "NAN", f17(apply(pp, 1, max)))))
cat("# largest per-unit range of p over the 8 seeds:", max(apply(pp, 1, function(v) diff(range(v))), na.rm = TRUE),
    " at counts", paste(a[which.max(apply(pp, 1, function(v) diff(range(v))))], b[which.max(apply(pp, 1, function(v) diff(range(v))))]), "\n")
# the NA defect (dossier sppt.md section 3): NA counts become 0 silently
dn <- d; dn$b[5] <- NA
rn <- suppressWarnings(sppt(dn, count_col = c("a", "b"), B = 999, seed = 1, arrangement = FALSE,
                            create_maps = FALSE, verbose = FALSE))
cat("# unit 5 with an NA test count: difference", rn$difference[5], " (base share", rn$base_prop[5],
    "): treated as 0 events, no warning\n")
set.seed(9); big <- st_sf(group = 1:10000, a = rpois(1e4, 5), b = rpois(1e4, 5),
                         geometry = st_make_grid(st_sfc(st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0))))), n = c(100, 100)))
t <- replicate(3, system.time(suppressWarnings(sppt(big, count_col = c("a", "b"), B = 999, seed = 1,
  arrangement = FALSE, create_maps = FALSE, verbose = FALSE)))[["elapsed"]])
cat(sprintf("TIME sppt 10,000 units B = 999: median %.3f s (k=3)\n", median(t)))
