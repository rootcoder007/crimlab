# Speed of the spatial functions against their counterparts, on fixed
# seeded point sets. Values must agree (stated tolerance) before a time is
# recorded. Writes inst/bench/sp.csv; nothing in it is typed by hand.
#   Rscript dev/bench-sp.R        (from the package root, crimlab installed)
# The table of record is taken on l14 (x86-64), one process, after the
# 1-minute load average has fallen below 1.0. "faster" compares with
# third-party packages only (base R, stats, spatstat, spdep, geosphere);
# none of the maintainer's own packages is a counterpart here.
suppressMessages({
  library(crimlab)
  library(spdep)
  library(spatstat.explore)
})
load <- function() as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
load0 <- load()
med <- function(f, k = 3) {
  median(vapply(seq_len(k), function(i) system.time(f())[["elapsed"]], 0))
}
pts <- function(n) {
  set.seed(n)
  cbind(x = runif(n), y = runif(n))
}
rows <- list()
add <- function(fn, ref, n, ok, tc, tr, note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    fn = fn, ref = ref, n = n, agree = ok, crimlab = tc, other = tr,
    ratio = if (ok) tr / tc else NA_real_,
    faster = if (ok) tr > tc else NA, note = note)
}
for (n in c(1e3, 1e4, 1e5)) {
  p <- pts(n)
  x <- p[, 1]
  y <- p[, 2]
  ok <- isTRUE(all.equal(unname(cgstat(p)[c("mx", "sx")]),
                         c(mean(x), sd(x)), tolerance = 1e-12))
  add("cgstat", "base mean, sd", n, ok, med(function() cgstat(p)),
      med(function() c(mean(x), mean(y), sd(x), sd(y))))
  e <- sdell(p)
  ok <- isTRUE(all.equal(unname(e[["a"]]),
                         sqrt(eigen(cov(p))$values[1]), tolerance = 1e-10))
  add("sdell", "base eigen(cov())", n, ok, med(function() sdell(p)),
      med(function() eigen(cov(p), symmetric = TRUE)))
  ok <- setequal(cvhull(p)$idx, grDevices::chull(p))
  add("cvhull", "grDevices chull", n, ok, med(function() cvhull(p)),
      med(function() grDevices::chull(p)))
  X <- spatstat.geom::ppp(x, y, c(0, 1), c(0, 1))
  ok <- isTRUE(all.equal(nnidx(p, win = c(0, 1, 0, 1), dist = TRUE)$dist[, 1],
                         spatstat.geom::nndist(X), tolerance = 1e-12))
  add("nnidx", "spatstat nndist", n, ok,
      med(function() nnidx(p, win = c(0, 1, 0, 1))),
      med(function() mean(spatstat.geom::nndist(X))))
  r <- seq(0, 2 / sqrt(n), length.out = 20)
  ok <- isTRUE(all.equal(ripk(p, win = c(0, 1, 0, 1), r = r)$k,
                         Kest(X, r = r, correction = "isotropic")$iso,
                         tolerance = 1e-10))
  add("ripk", "spatstat Kest iso", n, ok,
      med(function() ripk(p, win = c(0, 1, 0, 1), r = r)),
      med(function() Kest(X, r = r, correction = "isotropic")))
  s <- p[1:10, ]
  km <- kmpart(p, k = 10, init = s)
  b <- kmeans(p, s, iter.max = 100, algorithm = "Lloyd")
  ok <- identical(km$cl, unname(b$cluster))
  add("kmpart", "stats kmeans Lloyd", n, ok,
      med(function() kmpart(p, k = 10, init = s)),
      med(function() suppressWarnings(kmeans(p, s, iter.max = 100,
                                             algorithm = "Lloyd"))))
  ll <- cbind(lon = x * 360 - 180, lat = y * 180 - 90)
  ok <- isTRUE(all.equal(hvdist(ll, c(0, 0)), geosphere::distVincentySphere(
    ll, c(0, 0), r = 6371008.8), tolerance = 1e-12))
  add("hvdist", "geosphere distVincentySphere", n, ok,
      med(function() hvdist(ll, c(0, 0))),
      med(function() geosphere::distVincentySphere(ll, c(0, 0),
                                                   r = 6371008.8)))
  w <- spwt(p, k = 6)
  kn <- knn2nb(knearneigh(p, k = 6))
  lw <- nb2listw(kn)
  v <- x + y + rnorm(n)
  ok <- all(mapply(setequal, split(w$j, w$i), unclass(kn)))
  add("spwt", "spdep knearneigh + knn2nb + nb2listw", n, ok,
      med(function() spwt(p, k = 6)),
      med(function() nb2listw(knn2nb(knearneigh(p, k = 6)))))
  ok <- isTRUE(all.equal(gearc(v, w)[["c"]],
                         unname(geary.test(v, lw)$estimate[1]), 1e-10))
  add("gearc", "spdep geary.test", n, ok, med(function() gearc(v, w)),
      med(function() geary.test(v, lw)))
  pos <- exp(v)
  ok <- isTRUE(all.equal(gorg(pos, w)[["g"]],
                         unname(globalG.test(pos, lw)$estimate[1]), 1e-10))
  add("gorg", "spdep globalG.test", n, ok, med(function() gorg(pos, w)),
      med(function() globalG.test(pos, lw)))
  ok <- isTRUE(all.equal(lmorn(v, w)$ii,
                         unname(localmoran(v, lw, conditional = FALSE)[, 1]),
                         tolerance = 1e-10))
  add("lmorn", "spdep localmoran", n, ok, med(function() lmorn(v, w)),
      med(function() localmoran(v, lw, conditional = FALSE)))
  ok <- isTRUE(all.equal(locgi(v, w, star = FALSE),
                         as.vector(localG(v, lw)), tolerance = 1e-10))
  add("locgi", "spdep localG", n, ok, med(function() locgi(v, w, FALSE)),
      med(function() localG(v, lw)))
  if (n <= 1e4) {
    d <- as.matrix(dist(p))
    ok <- identical(fzmode(p, r = 0.02), as.integer(rowSums(d <= 0.02)))
    add("fzmode", "base dist matrix", n, ok, med(function() fzmode(p, r = 0.02)),
        med(function() rowSums(as.matrix(dist(p)) <= 0.02)))
  } else {
    add("fzmode", "base dist matrix", n, FALSE, NA, NA,
        "not run: the n x n distance matrix needs 80 GB at 1e5")
  }
}
out <- do.call(rbind, rows)
cpu <- sub(".*: ", "", grep("^model name", readLines("/proc/cpuinfo"),
                            value = TRUE)[1])
out$mach <- Sys.info()[["nodename"]]
out$cpu <- if (is.na(cpu)) Sys.info()[["machine"]] else cpu
out$rver <- R.version.string
out$load0 <- load0
out$load1 <- load()
out$date <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
# valid only when the machine was quiet: 1-minute load below 1.0 at the
# start (the caller waits for it) and below 1.5 at the end
out$valid <- load0 < 1 & out$load1 < 1.5
dir.create("inst/bench", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "inst/bench/sp.csv", row.names = FALSE)
print(out)
