# Reference values of tests/testthat/test-spw.R: spdep's neighbours,
# Geary's C, global and local G and local Moran, and circular's
# directional statistics, on the inputs the tests build. From the package
# root, with spdep and circular in the library:
#   Rscript dev/ref/spw.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("spdep", "circular")
source(file.path("dev", "ref", "put.R"))

# "weights agree with spdep's neighbours"
set.seed(31)
p <- cbind(x = runif(40), y = runif(40))
kn <- spdep::knn2nb(spdep::knearneigh(p, k = 3))
dn <- spdep::dnearneigh(p, 0, 0.2)
nbr <- list(
  inp = sum(p),
  kn = lapply(unclass(kn), sort),
  dn = lapply(1:40, function(i) setdiff(dn[[i]], 0L))
)

# "Geary's C, global and local G, local Moran agree with spdep"
set.seed(32)
p <- cbind(x = runif(50), y = runif(50))
v <- p[, 1] + p[, 2] + rnorm(50, sd = 0.3)
pos <- exp(2 * p[, 1]) + rexp(50)
lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(p, k = 4)))
r <- spdep::geary.test(v, lw)
n <- spdep::geary.test(v, lw, randomisation = FALSE)
dn <- spdep::dnearneigh(p, 0, 0.3)
lb <- spdep::nb2listw(dn, style = "B")
s <- spdep::globalG.test(pos, lb, zero.policy = TRUE)
sl <- spdep::localmoran(v, lw, conditional = FALSE)
ls <- spdep::localG(pos, spdep::nb2listw(spdep::include.self(dn),
                                         style = "B"))
l0 <- spdep::localG(pos, lb)
# the calls that made the weights would be evaluated when read back
attr(lw, "call") <- NULL
attr(lw$neighbours, "call") <- NULL
stat <- list(
  inp = c(sum(p), sum(v), sum(pos)),
  lw = lw,
  gc = unname(r$estimate[1]),
  gvr = unname(r$estimate[3]),
  gzr = unname(r$statistic),
  gpr = r$p.value,
  gvn = unname(n$estimate[3]),
  card = spdep::card(dn),
  gg = unname(s$estimate[1]),
  ge = unname(s$estimate[2]),
  gv = unname(s$estimate[3]),
  gz = as.vector(s$statistic),
  lmi = unname(sl[, 1]),
  lme = unname(sl[, 2]),
  lmv = unname(sl[, 3]),
  lmp = unname(sl[, 5]),
  gs = as.vector(ls),
  g0 = as.vector(l0)
)

# "directional means agree with circular"
a <- c(0.1, 0.4, 6.1, 1.2, 5.9)
cc <- circular::circular(a)
dir <- list(
  theta = as.numeric(circular::mean.circular(cc)) %% (2 * pi),
  r = circular::rho.circular(cc),
  var = circular::var.circular(cc),
  sd = circular::sd.circular(cc)
)

put(list(nbr = nbr, stat = stat, dir = dir), "spw", need)
