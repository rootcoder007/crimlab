# Speed of the spatial functions against their counterparts, on fixed
# seeded point sets. Values must agree (stated tolerance) before a time is
# recorded. Writes inst/bench/sp.csv; nothing in it is typed by hand.
#   Rscript dev/bench-sp.R        (from the package root, crimlab installed)
# The table of record is taken on one x86-64 Linux machine, one process, after the
# 1-minute load average has fallen below 1.0. "faster" compares with
# third-party packages only (base R, stats, spatstat, spdep, geosphere);
# none of the maintainer's own packages is a counterpart here.
# the 1-minute load average, from uptime (macOS and Linux)
load <- function() {
  v <- sub(".*load averages?: *", "", system("uptime", intern = TRUE))
  as.numeric(sub(",", ".", regmatches(v, regexpr("[0-9]+[.,][0-9]+", v))))
}
load0 <- load()
cat("crimlab from", find.package("crimlab"), "\n")
for (p in c("crimlab", "spdep", "spatstat.geom", "spatstat.explore",
            "geosphere")) {
  cat(p, format(packageVersion(p)), "\n")
}
suppressMessages({
  library(crimlab)
  library(spdep)
  library(spatstat.explore)
})
# seconds per call: each of k timings repeats f for at least 0.2 s
med <- function(f, k = 3) {
  median(vapply(seq_len(k), function(i) {
    r <- 0
    t0 <- proc.time()[["elapsed"]]
    repeat {
      f()
      r <- r + 1
      t <- proc.time()[["elapsed"]] - t0
      if (t >= 0.2) break
    }
    t / r
  }, 0))
}
pts <- function(n) {
  set.seed(n)
  cbind(x = runif(n), y = runif(n))
}
rows <- list()
# time one case (crimlab's call fc against the other's fr), recording the
# load before and after; a case whose load rose by more than 1.5 while it
# ran is timed again (up to three times)
add <- function(fn, ref, n, ok, fc, fr, note = "") {
  tc <- tr <- NA_real_
  l0 <- l1 <- load()
  tries <- 0L
  if (is.function(fc)) {
    repeat {
      tries <- tries + 1L
      l0 <- load()
      tc <- med(fc)
      tr <- med(fr)
      l1 <- load()
      if (l1 - l0 <= 1.5 || tries >= 3L) break
    }
  }
  rows[[length(rows) + 1L]] <<- data.frame(
    fn = fn, ref = ref, n = n, agree = ok, crimlab = tc, other = tr,
    ratio = if (ok) tr / tc else NA_real_,
    faster = if (ok) tr > tc else NA, note = note, lstart = l0, lend = l1,
    tries = tries)
}
for (n in c(1e3, 1e4, 1e5)) {
  p <- pts(n)
  x <- p[, 1]
  y <- p[, 2]
  ok <- isTRUE(all.equal(unname(cgstat(p)[c("mx", "sx")]),
                         c(mean(x), sd(x)), tolerance = 1e-12))
  add("cgstat", "base mean, sd", n, ok, (function() cgstat(p)),
      (function() c(mean(x), mean(y), sd(x), sd(y))))
  e <- sdell(p)
  ok <- isTRUE(all.equal(unname(e[["a"]]),
                         sqrt(eigen(cov(p))$values[1]), tolerance = 1e-10))
  add("sdell", "base eigen(cov())", n, ok, (function() sdell(p)),
      (function() eigen(cov(p), symmetric = TRUE)))
  ok <- setequal(cvhull(p)$idx, grDevices::chull(p))
  add("cvhull", "grDevices chull", n, ok, (function() cvhull(p)),
      (function() grDevices::chull(p)))
  X <- spatstat.geom::ppp(x, y, c(0, 1), c(0, 1))
  ok <- isTRUE(all.equal(nnidx(p, win = c(0, 1, 0, 1), dist = TRUE)$dist[, 1],
                         spatstat.geom::nndist(X), tolerance = 1e-12))
  add("nnidx", "spatstat nndist", n, ok,
      (function() nnidx(p, win = c(0, 1, 0, 1))),
      (function() mean(spatstat.geom::nndist(X))))
  r <- seq(0, 2 / sqrt(n), length.out = 20)
  ok <- isTRUE(all.equal(ripk(p, win = c(0, 1, 0, 1), r = r)$k,
                         Kest(X, r = r, correction = "isotropic")$iso,
                         tolerance = 1e-10))
  add("ripk", "spatstat Kest iso", n, ok,
      (function() ripk(p, win = c(0, 1, 0, 1), r = r)),
      (function() Kest(X, r = r, correction = "isotropic")))
  s <- p[1:10, ]
  km <- kmpart(p, k = 10, init = s)
  b <- kmeans(p, s, iter.max = 100, algorithm = "Lloyd")
  ok <- identical(km$cl, unname(b$cluster))
  add("kmpart", "stats kmeans Lloyd", n, ok,
      (function() kmpart(p, k = 10, init = s)),
      (function() suppressWarnings(kmeans(p, s, iter.max = 100,
                                             algorithm = "Lloyd"))))
  ll <- cbind(lon = x * 360 - 180, lat = y * 180 - 90)
  ok <- isTRUE(all.equal(hvdist(ll, c(0, 0)), geosphere::distVincentySphere(
    ll, c(0, 0), r = 6371008.8), tolerance = 1e-12))
  add("hvdist", "geosphere distVincentySphere", n, ok,
      (function() hvdist(ll, c(0, 0))),
      (function() geosphere::distVincentySphere(ll, c(0, 0),
                                                   r = 6371008.8)))
  w <- spwt(p, k = 6)
  kn <- knn2nb(knearneigh(p, k = 6))
  lw <- nb2listw(kn)
  v <- x + y + rnorm(n)
  ok <- all(mapply(setequal, split(w$j, w$i), unclass(kn)))
  add("spwt", "spdep knearneigh + knn2nb + nb2listw", n, ok,
      (function() spwt(p, k = 6)),
      (function() nb2listw(knn2nb(knearneigh(p, k = 6)))))
  ok <- isTRUE(all.equal(gearc(v, w)[["c"]],
                         unname(geary.test(v, lw)$estimate[1]), 1e-10))
  add("gearc", "spdep geary.test", n, ok, (function() gearc(v, w)),
      (function() geary.test(v, lw)))
  pos <- exp(v)
  ok <- isTRUE(all.equal(gorg(pos, w)[["g"]],
                         unname(globalG.test(pos, lw)$estimate[1]), 1e-10))
  add("gorg", "spdep globalG.test", n, ok, (function() gorg(pos, w)),
      (function() globalG.test(pos, lw)))
  ok <- isTRUE(all.equal(lmorn(v, w)$ii,
                         unname(localmoran(v, lw, conditional = FALSE)[, 1]),
                         tolerance = 1e-10))
  add("lmorn", "spdep localmoran", n, ok, (function() lmorn(v, w)),
      (function() localmoran(v, lw, conditional = FALSE)))
  ok <- isTRUE(all.equal(locgi(v, w, star = FALSE),
                         as.vector(localG(v, lw)), tolerance = 1e-10))
  add("locgi", "spdep localG", n, ok, (function() locgi(v, w, FALSE)),
      (function() localG(v, lw)))
  if (n <= 1e4) {
    d <- as.matrix(dist(p))
    ok <- identical(fzmode(p, r = 0.02), as.integer(rowSums(d <= 0.02)))
    add("fzmode", "base dist matrix", n, ok, (function() fzmode(p, r = 0.02)),
        (function() rowSums(as.matrix(dist(p)) <= 0.02)))
  } else {
    add("fzmode", "base dist matrix", n, FALSE, NA, NA,
        "not run: the n x n distance matrix needs 80 GB at 1e5")
  }
}
out <- do.call(rbind, rows)
# why a function is not faster, where it is not
why <- c(
  cgstat = "two passes with compensated sums, weights and geometric and harmonic means, against two plain sums",
  sdell = "compensated two-pass moments and the eigen-decomposition, against base's single cov() and eigen()",
  cvhull = "an exact orientation predicate (expansion arithmetic), so near-collinear points never turn wrong; chull uses a floating-point test",
  nnidx = "also computes the k-order expectations and the edge-corrected index",
  kmpart = "compensated sums for the centres on every pass and each group's sums of squares and of distances at the end; stats::kmeans' Lloyd is a plain Fortran loop")
slow <- !is.na(out$faster) & !out$faster
out$note[slow] <- paste0("not faster: ", why[out$fn[slow]])
cpu <- if (file.exists("/proc/cpuinfo")) {
  sub(".*: ", "", grep("^model name", readLines("/proc/cpuinfo"),
                       value = TRUE)[1])
} else {
  system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
}
out$mach <- Sys.info()[["nodename"]]
out$cpu <- if (is.na(cpu)) Sys.info()[["machine"]] else cpu
out$rver <- R.version.string
out$load0 <- load0
out$load1 <- load()
out$cores <- parallel::detectCores()
out$ram <- if (Sys.info()[["sysname"]] == "Darwin") {
  as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 2^30
} else NA_real_
out$date <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
# a case counts when the load did not rise by more than 1.5 while it ran
# (on the Mac of record the sessions that drive it keep the load near 2)
out$valid <- out$lend - out$lstart <= 1.5
dir.create("inst/bench", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "inst/bench/sp.csv", row.names = FALSE)
print(out)
