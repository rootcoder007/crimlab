# Reference values for src/trj.cpp (clab_zipll,
# clab_zipem). Field package: crimCV 1.0.0 (CRAN), R 4.5.0.
source("dev/doss/pre.R")
suppressMessages(library(crimCV))
cat("crimCV", format(packageVersion("crimCV")), "\n")
q <- function(e) { capture.output(r <- e); r }

# simulated two-group ZIP trajectories: 100 units x 6 periods
set.seed(11)
n <- 100; nt <- 6; oyr <- 0:(nt - 1) / (nt - 1)
X <- splines::bs(oyr, degree = 2, intercept = TRUE)   # 6 x 3
Z <- splines::bs(oyr, degree = 1, intercept = TRUE)   # 6 x 2
g <- rep(1:2, c(40, 60))
lam <- cbind(exp(0.2 + 1.5 * oyr), exp(1.2 - 1.0 * oyr))
rho <- c(0.25, 0.10)
Y <- t(sapply(1:n, function(i) ifelse(runif(nt) < rho[g[i]], 0, rpois(nt, lam[, g[i]]))))
Y[3, 4] <- -1                                         # one missing cell (crimCV's coding)
cat("\n// Y (100 x 6, row-major as listed), X = bs(oyr, 2, intercept), Z = bs(oyr, 1, intercept)\n")
cat(sprintf("const int Y[600] = {%s};\n", paste(t(Y), collapse = ", ")))
carr("Xb", X); carr("Zb", Z)

set.seed(3)
fit <- q(crimCV(Y, 2, dpolyp = 2, dpolyl = 1, model = "ZIP", init = 10))
cat("// crimCV(Y, 2, dpolyp = 2, dpolyl = 1, model = \"ZIP\", init = 10), set.seed(3)\n")
carr("cbeta", fit$beta); carr("cgamma", fit$gamma); carr("cprob", fit$prob)
cat(sprintf("const double cllike = %s;\n", f17(fit$llike)))
carr("cgwt", fit$gwt[1:5, ])

# Poisson (no zero inflation), one group: base R glm on the same cells
ok <- c(t(Y)) >= 0
yy <- c(t(Y))[ok]; XX <- X[rep(1:nt, n), ][ok, ]
gl <- glm(yy ~ XX - 1, family = poisson)
cat("// glm(y ~ X - 1, poisson) over the observed cells\n")
carr("pbeta", coef(gl)); cat(sprintf("const double pll = %s;\n", f17(as.numeric(logLik(gl)))))

# timing: ZIPt ng = 2, n = 1,000 x 10, init = 20 (default),
# and the ZIP model with the same data
set.seed(1); n2 <- 1000; nt2 <- 10; o2 <- 0:(nt2 - 1) / (nt2 - 1); g2 <- rbinom(n2, 1, 0.4) + 1
l2 <- cbind(exp(0.2 + 1.5 * o2), exp(1.2 - 1.0 * o2))
Y2 <- t(sapply(1:n2, function(i) ifelse(runif(nt2) < 0.15, 0, rpois(nt2, l2[, g2[i]]))))
t1 <- system.time(q(crimCV(Y2, 2, model = "ZIPt")))[["elapsed"]]
t2 <- system.time(f2 <- q(crimCV(Y2, 2, dpolyp = 3, dpolyl = 3, model = "ZIP")))[["elapsed"]]
cat(sprintf("\nTIME crimCV ZIPt ng=2 n=1000x10 init=20: %.3f s (1 run)\n", t1))
cat(sprintf("TIME crimCV ZIP ng=2 n=1000x10 dpolyp=3 dpolyl=3 init=20: %.3f s (1 run), llike %.6f\n", t2, f2$llike))
