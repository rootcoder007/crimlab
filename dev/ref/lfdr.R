# Reference values of tests/testthat/test-ibench.R: the local false
# discovery rate of locfdr::locfdr() (its defaults) on the z statistics of
# ibench() for the officers the test builds, the z statistics computed
# here with glm() and the propensity weights. From the package root:
#   Rscript dev/ref/lfdr.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- "locfdr"
source(file.path("dev", "ref", "put.R"))

set.seed(7)
no <- 150
off <- rep(sprintf("o%03d", seq_len(no)), each = 40)
night <- rbinom(no * 40, 1, rep(runif(no, 0.2, 0.8), each = 40))
bad <- off %in% sprintf("o%03d", 1:8)
y <- rbinom(no * 40, 1, plogis(-1 + night + 1.2 * bad))
z <- vapply(sort(unique(off)), function(o) {
  d <- as.numeric(off == o)
  e <- fitted(glm(d ~ night, family = binomial,
                  control = glm.control(epsilon = 1e-14, maxit = 100)))
  w <- (e / (1 - e))[d == 0]
  b <- sum(w * y[d == 0]) / sum(w)
  ess <- sum(w)^2 / sum(w^2)
  rt <- mean(y[d == 1])
  (rt - b) / sqrt(rt * (1 - rt) / sum(d) + b * (1 - b) / ess)
}, 0)
f <- locfdr::locfdr(unname(z), plot = 0)
put(list(z = unname(z), fdr = f$fdr, null = unname(f$fp0["mlest", 1:3])),
    "lfdr", need)
