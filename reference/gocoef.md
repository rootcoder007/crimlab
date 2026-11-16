# The slopes of a generalised ordered logit by threshold

The coefficients of
[`gologt()`](https://rootcoder007.github.io/crimlab/reference/gologt.md)
laid out as a table: one row per threshold of the outcome, one column
per term, so that a covariate's slopes can be read across the cuts of
the scale (a parallel covariate repeats its one slope), with their
standard errors and, for the free covariates, the Wald test that the
slopes are equal.

## Usage

``` r
gocoef(object, v = NULL)
```

## Arguments

- object:

  A result of
  [`gologt()`](https://rootcoder007.github.io/crimlab/reference/gologt.md).

- v:

  `NULL` for every term, or the names of terms: `"(Intercept)"` and the
  covariates.

## Value

A list: `est` and `se`, matrices with one row per threshold (named as
`"low|med"`) and one column per term of `v`; `wald`, the rows of the
fit's Wald tests for the covariates of `v`.

## See also

[`gologt()`](https://rootcoder007.github.io/crimlab/reference/gologt.md).

## Examples

``` r
set.seed(8)
n <- 300
d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
u <- runif(n)
d$risk <- 1 + (u < plogis(0.5 + 0.8 * d$race)) +
  (u < plogis(-1 + 0.1 * d$race - 0.3 * d$age))
gocoef(gologt(d, "risk", c("race", "age")), "race")
#> $est
#>            race
#> 1|2  0.69848026
#> 2|3 -0.03629893
#> 
#> $se
#>          race
#> 1|2 0.2765599
#> 2|3 0.2509516
#> 
#> $wald
#>   term     stat df          p
#> 1 race 6.354064  1 0.01171127
#> 
```
