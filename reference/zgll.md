# Group-based trajectory model for counts

The finite mixture of Nagin (1999) for counts over time, zero-inflated
Poisson within groups (Lambert 1992; Nielsen et al. 2014): unit `i`
belongs to group `k` with probability `prob[k]`; given `k`, its counts
are independent Poisson with log mean `off + x_t' beta_k` and, when `z`
is given, zero with extra probability `plogis(z_t' gamma_k)`. `zgll()`
evaluates the log-likelihood and the posterior memberships; `zgem()`
fits by a generalised EM from the values given (the log-likelihood never
decreases). One start per call: the model has local maxima, so several
starts are the caller's loop, and the number of groups is the caller's
choice (compare the log-likelihoods with a penalty such as BIC).

## Usage

``` r
zgll(y, x, z = NULL, par, prob, off = NULL)

zgem(y, x, z = NULL, par, prob, off = NULL, maxit = 500, tol = 1e-08)
```

## Arguments

- y:

  Counts: a matrix with one row per unit and one column per period; `NA`
  or a negative value is a missing period.

- x:

  The design of the Poisson part: a matrix with one row per period
  (columns 1, t, t^2, say).

- z:

  `NULL` for a Poisson model, or the design of the zero-inflation part,
  one row per period.

- par:

  The coefficients: a matrix with one row per column of `x` (then of
  `z`) and one column per group.

- prob:

  The group shares: positive numbers summing to 1.

- off:

  `NULL`, or log exposures with the dimensions of `y`.

- maxit:

  The most iterations: a whole number from 1.

- tol:

  Stop when an iteration gains at most `tol (1 + |ll|)`.

## Value

`zgll()`: a list of `ll` and `post` (a matrix of posterior memberships,
one row per unit). `zgem()`: a list of the fitted `par`, `prob`, `post`,
`ll` and `iter` (iterations made).

## References

Nagin, D. S. (1999). Analyzing developmental trajectories: a
semiparametric, group-based approach. *Psychological Methods*, 4(2),
139-157.
[doi:10.1037/1082-989X.4.2.139](https://doi.org/10.1037/1082-989X.4.2.139)

Lambert, D. (1992). Zero-inflated Poisson regression, with an
application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
[doi:10.2307/1269547](https://doi.org/10.2307/1269547)

## Examples

``` r
set.seed(1)
x <- cbind(1, (1:6 - 3.5) / 2.5)
b <- cbind(c(0.5, -1), c(1.5, 0.5))
y <- t(sapply(rep(1:2, each = 50),
              function(k) rpois(6, exp(x %*% b[, k]))))
f <- zgem(y, x, par = b + 0.2, prob = c(0.5, 0.5))
f$par
#>            [,1]     [,2]
#> [1,]  0.4069311 1.530861
#> [2,] -1.1301323 0.481301
```
