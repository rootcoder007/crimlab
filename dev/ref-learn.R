# Reference values for tests/testthat/test-learn.R, from glmnet and glm():
#   Rscript dev/ref-learn.R      (from the package root; needs glmnet)
# Writes tests/testthat/ref-learn.rds: the data, the folds, and what
# glmnet (and cv.glmnet) and glm() return on them. Both sides run with the
# convergence threshold 1e-20 (glmnet's control list), so they agree to far below the tolerance the
# test states (1e-8 on coefficients).
stopifnot(requireNamespace("glmnet"))
set.seed(20261112)
n <- 150
p <- 8
z <- matrix(rnorm(n * p), n)
x <- z %*% chol(0.4^abs(outer(seq_len(p), seq_len(p), "-")))
x <- sweep(x, 2, c(1, 3, 0.5, 2, 1, 10, 1, 0.2), "*")
x <- sweep(x, 2, c(0, 5, -2, 0, 1, 100, 0, 3), "+")
colnames(x) <- paste0("v", seq_len(p))
b <- c(1, -0.5, 2, 0, 0, 0.05, 0, 1)
eta <- drop(scale(x) %*% b)
y <- 2 + eta + rnorm(n)
yb <- rbinom(n, 1, plogis(0.3 + eta))
fold <- sample(rep(1:5, length.out = n))
ctl <- list(thresh = 1e-20, maxit = 1e7)
one <- function(fam, yy, alpha) {
  # glmnet's own path; then the fit at those lambdas given explicitly:
  # on its own path glmnet reports zero coefficients at the first lambda
  # even for alpha = 0, where the ridge solution is not zero
  lam <- glmnet::glmnet(x, yy, family = fam, alpha = alpha, nlambda = 15,
                        lambda.min.ratio = 0.02, control = ctl)$lambda
  f <- glmnet::glmnet(x, yy, family = fam, alpha = alpha, lambda = lam,
                      control = ctl)
  cv <- glmnet::cv.glmnet(x, yy, family = fam, alpha = alpha,
                          lambda = f$lambda, foldid = fold, control = ctl)
  list(fam = fam, alpha = alpha, lambda = f$lambda, a0 = unname(f$a0),
       beta = unname(as.matrix(f$beta)), cvm = cv$cvm, cvsd = cv$cvsd,
       lmin = cv$lambda.min, l1se = cv$lambda.1se)
}
fits <- list(
  one("gaussian", y, 1), one("gaussian", y, 0.5), one("gaussian", y, 0),
  one("binomial", yb, 1), one("binomial", yb, 0.3)
)
g <- glm(yb ~ x, family = binomial(),
         control = glm.control(epsilon = 1e-14, maxit = 100))
ref <- list(x = x, y = y, yb = yb, fold = fold, fits = fits,
            glm = list(coef = unname(coef(g)),
                       se = unname(sqrt(diag(vcov(g)))),
                       dev = deviance(g)),
            made = list(R = R.version.string,
                        glmnet = as.character(packageVersion("glmnet"))))
saveRDS(ref, "tests/testthat/ref-learn.rds", version = 2)
