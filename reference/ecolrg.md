# Ecological regression with expected counts

Regression of area counts `y` on covariates with the offset `log(e)`,
`e` the expected counts (from indirect standardisation, say), so that
the coefficients act on the standardised ratio: Poisson by iteratively
reweighted least squares; negative binomial, alternating Newton steps on
the score of `theta` with reweighted fits of the coefficients; or
zero-inflated Poisson (Lambert 1992) by EM, the zero part logistic on an
intercept and the columns of `z`.

## Usage

``` r
ecolrg(y, e, x = NULL, family = "poisson", z = NULL, tol = 1e-12, maxit = 500)
```

## Arguments

- y:

  The counts: numbers from 0.

- e:

  The expected counts: positive numbers, one per count.

- x:

  `NULL`, or the covariates: a numeric matrix (or vector) with one row
  per count. An intercept is added.

- family:

  `"poisson"`, `"negbin"` or `"zip"`.

- z:

  For `"zip"`: `NULL`, or the covariates of the zero part, as `x`.

- tol:

  The convergence tolerance of `theta` and of the zero-inflated
  log-likelihood: one number strictly between 0 and 1.

- maxit:

  The rounds: one whole number from 1.

## Value

A list of class `ecolrg`: `coef` (the intercept first), `theta`
(negative binomial), `zcoef` (zero-inflated), `fitted` (the means of the
count part) and `loglik`.

## References

Lambert, D. (1992). Zero-inflated Poisson regression, with an
application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
[doi:10.2307/1269547](https://doi.org/10.2307/1269547)

Lawson, A. B. (2018). *Bayesian Disease Mapping*, third edition. Chapman
and Hall/CRC, Boca Raton.
[doi:10.1201/9781351271769](https://doi.org/10.1201/9781351271769)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
# simulated counts and expected counts
set.seed(1)
x <- rnorm(30)
e <- runif(30, 3, 8)
y <- rpois(30, e * exp(0.2 + 0.4 * x))
ecolrg(y, e, x)
#> ecolrg: 30 counts, log-likelihood -67.3106
#> [1] 0.1404659 0.4851843
```
