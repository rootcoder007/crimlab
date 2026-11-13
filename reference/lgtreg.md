# Logistic regression

Maximum likelihood for a 0/1 outcome and a logit link, by Newton's
method on the design with centred columns, the method
[`glm()`](https://rdrr.io/r/stats/glm.html) uses with
`family = binomial()` when the columns are not collinear. Standard
errors come from the inverse Fisher information.

## Usage

``` r
lgtreg(x, y, tol = 1e-10, maxit = 50)
```

## Arguments

- x:

  The predictors: a numeric matrix with more rows than columns plus one,
  with no missing values; an intercept is always added.

- y:

  The outcome: 0 or 1, one per row of `x`, both present.

- tol:

  Newton steps stop when none moves a coefficient by more than
  `tol (1 + |coefficient|)`: one positive number.

- maxit:

  The most Newton steps: one whole number.

## Value

A list: `coef` and `se`, named numeric vectors with the intercept first,
and `dev`, the residual deviance.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models* (2nd
ed.). Chapman and Hall.
[doi:10.1007/978-1-4899-3242-6](https://doi.org/10.1007/978-1-4899-3242-6)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(300 * 2), 300)
y <- rbinom(300, 1, plogis(0.5 + x[, 1] - x[, 2]))
lgtreg(x, y)$coef
#> (Intercept)          x1          x2 
#>   0.4864763   1.0070780  -0.8877568 
```
