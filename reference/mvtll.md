# Group-based multivariate trajectories

The group-based trajectory model with several normal indicators per unit
(Nagin 1999; Magrini 2022): given group `k`, the `r` indicators of a
unit at period `t` are multivariate normal with mean `B_k' x_t` and
covariance `S_k`, independently over periods. `mvtll()` evaluates the
log-likelihood and posterior memberships; `mvtem()` fits by EM with
closed-form steps from the values given. A period with any missing
indicator is left out for that unit.

## Usage

``` r
mvtll(y, x, beta, sig, prob)

mvtem(y, x, beta, sig, prob, maxit = 500, tol = 1e-08)
```

## Arguments

- y:

  The indicators: an array of units by periods by indicators; `NA`
  allowed, not infinite values.

- x:

  The design: a matrix with one row per period.

- beta:

  The coefficients: an array of design columns by indicators by groups.

- sig:

  The covariances: an array of indicators by indicators by groups, each
  positive definite.

- prob:

  The group shares: positive numbers summing to 1.

- maxit:

  The most iterations: a whole number from 1.

- tol:

  Stop when an iteration gains at most `tol (1 + |ll|)`.

## Value

`mvtll()`: a list of `ll` and `post`. `mvtem()`: a list of the fitted
`beta`, `sig`, `prob`, `post`, `ll` and `iter`.

## References

Magrini, A. (2022). Assessment of agricultural sustainability in
European Union countries: a group-based multivariate trajectory
approach. *AStA Advances in Statistical Analysis*, 106.
[doi:10.1007/s10182-022-00437-9](https://doi.org/10.1007/s10182-022-00437-9)

## Examples

``` r
set.seed(1)
x <- cbind(1, 1:5)
y <- array(rnorm(40 * 5 * 2), c(40, 5, 2))
y[1:20, , 1] <- y[1:20, , 1] + 3
b <- array(c(0, 0, 0, 0, 3, 0, 0, 0), c(2, 2, 2))
s <- array(diag(2), c(2, 2, 2))
mvtem(y, x, b, s, prob = c(0.5, 0.5))$prob
#> [1] 0.5000024 0.4999976
```
