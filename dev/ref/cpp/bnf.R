# Reference values for src/bnf.cpp (clab_bdig,
# clab_bnft). Field packages: benford.analysis 0.1.5 and BenfordTests 1.2.0
# (CRAN), R 4.5.0.
source("dev/ref/cpp/pre.R")
suppressMessages({library(benford.analysis); library(BenfordTests)})
cat("benford.analysis", format(packageVersion("benford.analysis")),
    "BenfordTests", format(packageVersion("BenfordTests")), "\n")

# 1. deterministic first-digit sample (an integer rule), rebuilt in src/t_bnf.cpp
i <- 0:999
d1 <- 1 + ((i * 7 + (i * i) %% 13) %% 9)
x1 <- d1 + 0.5
ct <- chisq.benftest(x1)
cat("\n// BenfordTests::chisq.benftest(d + 0.5), d = 1 + ((i * 7 + (i^2 %% 13)) %% 9), i = 0..999\n")
cat(sprintf("const double ref_chi1[2] = {%s};  // statistic, p-value\n", f17(c(ct$statistic, ct$p.value))))
cat(sprintf("const double ref_cnt1[9] = {%s};\n", f17(tabulate(d1, 9))))
# a Benford-like sample (counts proportional to log10(1 + 1/d), rounded)
nb <- round(10000 * log10(1 + 1 / (1:9)))
xb <- rep(1:9, nb) + 0.5
cb <- chisq.benftest(xb)
cat("// nb = round(10000 log10(1 + 1/d)): chisq.benftest(rep(1:9, nb) + 0.5)\n")
cat(sprintf("const double ref_chib[2] = {%s};\n", f17(c(cb$statistic, cb$p.value))))

# 2. two digits, benford.analysis (sign = "both", so nothing is dropped)
d2 <- 10 + ((i * 37 + (i * i) %% 31) %% 90)
x2 <- d2 + 0.5
b2 <- benford(x2, number.of.digits = 2, sign = "both")
cat("\n// benford(d2 + 0.5, 2, sign = \"both\"), d2 = 10 + ((i * 37 + (i^2 %% 31)) %% 90)\n")
cat(sprintf("const double ref_chi2[3] = {%s};  // chisq statistic, df, p\n",
    f17(c(chisq(b2)$statistic, chisq(b2)$parameter, chisq(b2)$p.value))))
cat(sprintf("const double ref_mad2 = %s;  // MAD(b2)\n", f17(MAD(b2))))
t2 <- chisq.benftest(x2, digits = 2)
cat(sprintf("const double ref_chi2t[2] = {%s};  // BenfordTests chisq.benftest(x2, 2)\n",
    f17(c(t2$statistic, t2$p.value))))

# timing
set.seed(5); big <- rbenf(1e6)
t <- replicate(3, system.time(benford(big, number.of.digits = 2))[["elapsed"]])
cat(sprintf("\nTIME benford.analysis::benford n=1e6 2 digits: median %.3f s (k=3)\n", median(t)))
t <- replicate(3, system.time(chisq.benftest(big))[["elapsed"]])
cat(sprintf("TIME BenfordTests::chisq.benftest n=1e6: median %.3f s (k=3)\n", median(t)))
