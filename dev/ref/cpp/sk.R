# Reference values for src/sk.cpp (clab_skp). Field
# package: SK4FGA 0.1.1 (CRAN), R 4.5.0.
source("dev/ref/cpp/pre.R")
library(SK4FGA)
cat("SK4FGA", format(packageVersion("SK4FGA")), "\n")
# deterministic refractive indices, rebuilt in src/t_sk.cpp: three groups
# of 8 (means 1.5180, 1.5186, 1.5200), offsets ((j * 7) %% 9 - 4) * 1e-5
j <- 0:23
ri <- c(1.5180, 1.5186, 1.5200)[j %/% 8 + 1] + ((j * 7) %% 9 - 4) * 1e-5
cat("\npartition(ri)$groups:", partition(ri)$groups, "\n")
h <- 1.5183 + ((j[1:20] * 7) %% 9 - 4) * 1e-5
cat("partition(homogeneous 20)$groups:", partition(h)$groups, "\n")
cat("partition((ri - 1.5) * 1000)$groups:", partition((ri - 1.5) * 1000)$groups, "\n")
set.seed(2); x <- rnorm(5000, 1.518, 1e-4)
t <- replicate(3, system.time(partition(x))[["elapsed"]])
cat(sprintf("\nTIME partition n=5000: median %.3f s (k=3), %d groups\n", median(t),
            max(partition(x)$groups)))
