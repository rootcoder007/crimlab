# Count regression with a random intercept per cluster

Fits the log-linear model of a count with a normal random intercept for
each cluster, \\\log E\[y\_{ij} \| u_j\] = x\_{ij}'\beta + u_j\\, \\u_j
\sim N(0, \sigma^2)\\, with a Poisson or a negative binomial (variance
`mu + mu^2 / theta`) response, by maximum likelihood with the Laplace
approximation of each cluster's integral (Breslow and Clayton 1993), as
'glmmTMB' and 'lme4' fit it. It is the main model of Ruhela (2026):
placement volatility on alert complexity with a random intercept per
regional cluster.

## Usage

``` r
cntmix(y, x, grp, fam = c("nb", "poisson"), w = NULL, level = 0.95)
```

## Arguments

- y:

  The counts: whole numbers from 0.

- x:

  The covariates: a data frame (factors and text enter as indicators of
  each level but the first) or a numeric matrix, one row per count; an
  intercept is added.

- grp:

  The cluster of each count: labels, at least two clusters, one of them
  with two rows or more.

- fam:

  `"nb"` or `"poisson"`.

- w:

  `NULL` or a weight from 0 for each count.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `cntmix`, a list: `coef`, a data frame (`term`,
`est`, `se`, `z`, `p`, `irr`, `lo`, `hi`; `NA` for a term that cannot be
estimated); `sigma`, the standard deviation of the intercepts; `theta`
(`NA` for Poisson); `loglik`, `aic`; `fam`, the family reported; `n`,
`ngrp`; `conv`, the optimiser's code (0 when it converged); `level`;
`notes`, why a term has no estimate (empty when every term has one).

## Details

The algorithm: for given `beta`, `log sigma` (and `log theta`), each
cluster's mode of the random intercept is found by Newton's method on
its one-dimensional log-density, all clusters at once (`O(n)` per step,
quadratic convergence); the Laplace approximation of the marginal
log-likelihood is then a sum over clusters. The outer parameters are
maximised by BFGS from the Poisson fit without random effects, and their
standard errors are from the inverse of the numerically differentiated
Hessian at the maximum.

When the negative binomial's `theta` runs off to the Poisson limit
(above `1e6`), its likelihood is flat in that direction and the Hessian
is singular; the Poisson model is then reported, as it is the limit,
with `theta` noted (`fam` says which). Fewer than two clusters, or
clusters of one row each, are refused: the variance of the intercepts is
then not identified.

A coefficient that cannot be estimated is `NA` with its standard error,
test and interval, and `notes` says why, never a standard error of 0: a
column aliased with the others (a constant, or a combination of other
columns, found by the QR decomposition of the design before the fit,
which drops it as [`stats::lm()`](https://rdrr.io/r/stats/lm.html)
does); a 0/1 column all of whose counts on one side are 0 (separation:
its rate ratio is 0 or infinite, and so is the intercept's when the
zeros are where the column is 0); and an estimate whose variance from
the Hessian is not a positive number.

## References

Breslow, N. E. and Clayton, D. G. (1993). Approximate inference in
generalized linear mixed models. *Journal of the American Statistical
Association*, 88(421), 9-25.
[doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284)

Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg,
C. W., Nielsen, A., Skaug, H. J., Maechler, M. and Bolker, B. M. (2017).
glmmTMB balances speed and flexibility among packages for zero-inflated
generalized linear mixed modeling. *The R Journal*, 9(2), 378-400.
[doi:10.32614/RJ-2017-066](https://doi.org/10.32614/RJ-2017-066)

## Examples

``` r
set.seed(2)
g <- rep(1:12, each = 40)
d <- rbinom(480, 1, 0.3)
y <- rpois(480, exp(-1 + 0.4 * d + rnorm(12, 0, 0.5)[g]))
cntmix(y, data.frame(d = d), g, fam = "poisson")
#> cntmix: poisson counts with a random intercept per cluster
#> 480 counts in 12 clusters; sigma 0.5101; log-likelihood -459.3308
#>         term     est     se      z         p    irr     lo    hi
#>  (Intercept) -0.8465 0.1695 -4.994 5.919e-07 0.4289 0.3077 0.598
#>            d  0.4339 0.1215  3.573 3.535e-04 1.5433 1.2164 1.958
```
