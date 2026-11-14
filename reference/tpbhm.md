# Bayesian preservation model

A Bayesian linear model of a preservation score (or any outcome) on
burial conditions and measurements, with a normal prior on each
coefficient and the noise variance fixed at the least-squares residual
variance (an empirical Bayes plug-in), so that the posterior is normal
in closed form. An informative prior on a treatment (say, that quicklime
preserves) is set through `pm` and `psd`, and the data revise it. With
groups (cemeteries, sites) the mean residual of each group is shrunk
towards 0 by the method of moments: `tau2` is the variance between the
group means beyond what the noise explains.

## Usage

``` r
tpbhm(y, x, grp = NULL, pm = 0, psd = 10, level = 0.95)
```

## Arguments

- y:

  The outcome: finite numbers, at least two cases. A `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too (the `taph-*` presets): its elements fill the arguments
  by name.

- x:

  The predictors: a numeric matrix or data frame of finite numbers, one
  row per case (code a factor as indicator columns first).

- grp:

  `NULL`, or the group of each case: labels, at least two groups.

- pm, psd:

  The prior means and standard deviations of the intercept and the
  coefficients: one number each or one per term (`psd` positive).

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

A list: `coef`, a data frame with one row per term (`term`, `mean`,
`sd`, the interval `lo`, `hi`, `ppos` the posterior probability that the
coefficient is above 0, `prior` the prior's standard deviation);
`sigma`, the noise standard deviation; `tau2`; `groups`, `NULL` or a
data frame with one row per group (`group`, `n`, `raw` the mean
residual, `shrink` the factor, `pooled` the shrunk intercept); `fitted`;
`level`.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd ed. CRC Press.
[doi:10.1201/b16018](https://doi.org/10.1201/b16018)

## Examples

``` r
set.seed(1)
d <- data.frame(lime = rbinom(30, 1, 0.5), temp = rnorm(30, 12, 3))
y <- 0.3 + 0.2 * d$lime - 0.01 * d$temp + rnorm(30, sd = 0.05)
tpbhm(y, d, grp = rep(c("a", "b", "c"), 10), pm = c(0, 0.2, 0))$coef
#>          term        mean          sd          lo          hi      ppos prior
#> 1 (Intercept) 0.173071063 0.045858883  0.08318930 0.262952821 0.9999197    10
#> 2        lime 0.184805792 0.016328473  0.15280257 0.216809011 1.0000000    10
#> 3        temp 0.001683087 0.003552176 -0.00527905 0.008645224 0.6821854    10
```
