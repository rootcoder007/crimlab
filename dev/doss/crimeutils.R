# Reference values for src/rpt.cpp (clab_rptc).
# Field package: crimeutils 0.5.1 (CRAN). It does not install on zeus (gt
# needs V8), so its file R/ucr_constant_reporters.R is sourced from the
# CRAN tarball unpacked in /mnt/nvme/work/csls/src/crimeutils, R 4.5.0.
source("dev/doss/pre.R")
suppressMessages({library(dplyr); library(magrittr)})
e <- new.env()
sys.source("/mnt/nvme/work/csls/src/crimeutils/R/ucr_constant_reporters.R", e)
f <- e$ucr_constant_reporter_oris

# deterministic panel, rebuilt in src/t_rpt.cpp: agencies 0..39, years
# 2000..2009, months missing ((u * 7 + y * 3 + (u * y) %% 5) %% 13) for
# unit u and year index y, except that unit u reports nothing in year y
# when (u + 2 * y) %% 17 == 0
u <- rep(0:39, each = 10); y <- rep(0:9, 40)
miss <- (u * 7 + y * 3 + (u * y) %% 5) %% 13
keep <- (u + 2 * y) %% 17 != 0
p <- data.frame(ori = sprintf("u%02d", u), year = 2000 + y, number_of_months_missing = miss)[keep, ]
cat("\n// crimeutils 0.5.1 ucr_constant_reporter_oris(p, m): constant units\n")
for (m in c(9, 12)) {
  r <- sort(as.integer(sub("u", "", f(p, m))))
  cat(sprintf("// m = %d: {%s}\n", m, paste(r, collapse = ", ")))
}

set.seed(1); na <- 20000; ny <- 30
big <- data.frame(ori = rep(sprintf("a%05d", 1:na), each = ny), year = rep(1991:2020, na),
                  number_of_months_missing = rbinom(na * ny, 12, 0.02))
t <- replicate(3, system.time(f(big, 0))[["elapsed"]])
cat(sprintf("\nTIME ucr_constant_reporter_oris 20,000 agencies x 30 years: median %.3f s (k=3), %d constant\n",
    median(t), length(f(big, 0))))
