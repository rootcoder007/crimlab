# Reference values and defect reproducer for src/sk.cpp (clab_skp). Field
# package: SK4FGA 0.1.1 (CRAN), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/SK4FGA.R
library(SK4FGA)
cat("SK4FGA", format(packageVersion("SK4FGA")), "\n")
# deterministic refractive indices, rebuilt in src/t_sk.cpp: three groups
# of 8 (means 1.5180, 1.5186, 1.5200), offsets ((j * 7) %% 9 - 4) * 1e-5
j <- 0:23
ri <- c(1.5180, 1.5186, 1.5200)[j %/% 8 + 1] + ((j * 7) %% 9 - 4) * 1e-5
cat("\npartition(ri)$groups:", partition(ri)$groups, "\n")
h <- 1.5183 + ((j[1:20] * 7) %% 9 - 4) * 1e-5
cat("partition(homogeneous 20)$groups:", partition(h)$groups, "\n")
# defect 1 (dossier SK4FGA.md section 3): the variance 1.6e-9 is fixed
cat("partition((ri - 1.5) * 1000)$groups:", partition((ri - 1.5) * 1000)$groups, "\n")
# defect 2: NA dropped, labels no longer aligned with the input
g <- partition(c(ri[1:3], NA, ri[4:24]))$groups
cat("with an NA as 4th of 25 inputs:", length(g), "labels:", g, "\n")
# defect 3 (new): above k = 20 the chi-square threshold over-splits; size
# of the first split under the null, 20,000 samples each
set.seed(1); c0 <- pi / (2 * (pi - 2))
for (k in c(20, 21, 50, 200)) {
  lam <- replicate(20000, c0 * find_B0(sort(rnorm(k, 0, 4e-5)))$x / 1.6e-9)
  cat("k =", k, "size at alpha 0.05:", mean(lam > SK4FGA:::calculate_lambda_threshold(k, 0.05)), "\n")
}
set.seed(2); x <- rnorm(5000, 1.518, 1e-4)
t <- replicate(3, system.time(partition(x))[["elapsed"]])
cat(sprintf("\nTIME partition n=5000: median %.3f s (k=3), %d groups\n", median(t),
            max(partition(x)$groups)))
