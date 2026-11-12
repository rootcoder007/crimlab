# Reference values and defect reproducer for src/trj.cpp (clab_zipll,
# clab_zipem). Field package: crimCV 1.0.0 (CRAN), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/crimCV.R
suppressMessages(library(crimCV))
cat("crimCV", format(packageVersion("crimCV")), "\n")
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
carr <- function(nm, x) cat(sprintf("const double %s[%d] = {%s};\n", nm, length(x), f17(x)))
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

# defects (dossier crimCV.md section 3)
Yn <- Y; Yn[5, 2] <- NA
cat("\n# NA in one cell:", tryCatch({q(crimCV(Yn, 2, dpolyp = 2, dpolyl = 1, model = "ZIP", init = 2)); "ran"},
    error = conditionMessage), "\n")
cat("# n = 2 units:", tryCatch({q(crimCV(Y[1:2, ], 2, dpolyp = 2, dpolyl = 1, model = "ZIP", init = 2)); "ran"},
    error = conditionMessage), "\n")

# timing: the dossier's case, ZIPt ng = 2, n = 1,000 x 10, init = 20 (default),
# and the ZIP model with the same data
set.seed(1); n2 <- 1000; nt2 <- 10; o2 <- 0:(nt2 - 1) / (nt2 - 1); g2 <- rbinom(n2, 1, 0.4) + 1
l2 <- cbind(exp(0.2 + 1.5 * o2), exp(1.2 - 1.0 * o2))
Y2 <- t(sapply(1:n2, function(i) ifelse(runif(nt2) < 0.15, 0, rpois(nt2, l2[, g2[i]]))))
t1 <- system.time(q(crimCV(Y2, 2, model = "ZIPt")))[["elapsed"]]
t2 <- system.time(f2 <- q(crimCV(Y2, 2, dpolyp = 3, dpolyl = 3, model = "ZIP")))[["elapsed"]]
cat(sprintf("\nTIME crimCV ZIPt ng=2 n=1000x10 init=20: %.3f s (1 run)\n", t1))
cat(sprintf("TIME crimCV ZIP ng=2 n=1000x10 dpolyp=3 dpolyl=3 init=20: %.3f s (1 run), llike %.6f\n", t2, f2$llike))

# check (base R, independent of crimlab): the log-likelihood at the point
# clab_zipem reaches from crimCV's estimate (parameters rounded to 5
# decimals, from src/t_trj.cpp), against crimCV's reported maximum
llz <- function(beta, gamma, prob) {
  sum(sapply(1:nrow(Y), function(i) {
    o <- Y[i, ] >= 0
    log(sum(sapply(1:2, function(k) {
      lam <- exp(X %*% beta[, k]); rho <- plogis(Z %*% gamma[, k])
      prob[k] * prod((Y[i, o] == 0) * rho[o] + (1 - rho[o]) * dpois(Y[i, o], lam[o]))
    })))
  }))
}
cat("\n# base-R log-likelihood at crimCV's estimate:", llz(fit$beta, fit$gamma, fit$prob), "\n")
cat("# base-R log-likelihood at the EM point:",
    llz(matrix(c(0.17000, 1.18680, 1.60814, 1.24426, 0.63853, 0.25274), 3),
        matrix(c(-1.30619, -0.64623, -1.99047, -2.34619), 2), c(0.4201, 0.5799)), "\n")
