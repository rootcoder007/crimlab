# Reference values for tests/testthat/test-rf.R, from ranger:
#   Rscript dev/ref-rf.R      (from the package root; needs ranger)
# Known regression and probability functions, training and test rows;
# ranger's test error against the true function with its defaults
# (500 trees, mtry floor(sqrt(p)), node size 5, or 10 for probabilities),
# averaged over five seeds. Writes tests/testthat/ref-rf.rds.
stopifnot(requireNamespace("ranger"))
set.seed(20261112)
n <- 800
nt <- 1000
p <- 6
f <- function(x) 2 * sin(x[, 1]) + (x[, 2] > 0) - 0.5 * x[, 3]^2 + x[, 4]
x <- matrix(rnorm(n * p), n, dimnames = list(NULL, paste0("x", 1:p)))
xt <- matrix(rnorm(nt * p), nt, dimnames = list(NULL, paste0("x", 1:p)))
y <- f(x) + rnorm(n, sd = 0.5)
yb <- rbinom(n, 1, plogis(f(x)))
mse <- function(a, b) mean((a - b)^2)
reg <- mean(vapply(1:5, function(s) {
  r <- ranger::ranger(x = x, y = y, seed = s, num.threads = 1)
  mse(predict(r, xt)$predictions, f(xt))
}, 0))
prob <- mean(vapply(1:5, function(s) {
  r <- ranger::ranger(x = x, y = factor(yb), probability = TRUE, seed = s,
                      num.threads = 1)
  mse(predict(r, xt)$predictions[, "1"], plogis(f(xt)))
}, 0))
saveRDS(list(x = x, xt = xt, y = y, yb = yb, ft = f(xt),
             pt = plogis(f(xt)), reg = reg, prob = prob,
             made = list(R = R.version.string,
                         ranger = as.character(packageVersion("ranger")))),
        "tests/testthat/ref-rf.rds", version = 2)
