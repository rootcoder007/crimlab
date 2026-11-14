# Synthetic exposure of areas from a survey

A small-area synthetic estimate (Gonzalez and Hoza 1978): the prevalence
of a trait (serious mental illness, say, for the benchmark of a
disparity in police use of force) is modelled on survey microdata by
logistic regression
([`lgtreg()`](https://rootcoder007.github.io/crimlab/reference/lgtreg.md));
applied to each area's covariates it gives the area's rate, and times
its population the number of people with the trait, the exposure a rate
per person needs. The area covariates must be on the scale of the
survey's (shares of an area for indicators read as the probability of
the indicator, a plug-in that assumes the model holds within areas).

## Usage

``` r
synexp(y, x, ax, pop)
```

## Arguments

- y:

  The trait of each survey respondent: 0 or 1, both present.

- x:

  The respondents' covariates: a numeric matrix or data frame, one row
  per respondent, more rows than columns plus one.

- ax:

  The areas' covariates: a numeric matrix or data frame with the columns
  of `x`, one row per area.

- pop:

  The population of each area: numbers from 0.

## Value

A list: `coef`, the logistic fit (`term`, `est`, `se`); `areas`, a data
frame with one row per area (`rate`, `pop`, `expo`); `total`, the summed
exposure; `dev`, the fit's deviance.

## References

Gonzalez, M. E. and Hoza, C. (1978). Small-area estimation with
application to unemployment and housing estimates. *Journal of the
American Statistical Association*, 73(361), 7-15.
[doi:10.1080/01621459.1978.10479991](https://doi.org/10.1080/01621459.1978.10479991)

## Examples

``` r
set.seed(1)
s <- data.frame(age = runif(400, 18, 80), male = rbinom(400, 1, 0.5))
y <- rbinom(400, 1, plogis(-3 + 0.02 * s$age + 0.4 * s$male))
areas <- data.frame(age = c(35, 52), male = c(0.48, 0.51))
synexp(y, s, areas, pop = c(12000, 8000))$areas
#>         rate   pop     expo
#> 1 0.09284222 12000 1114.107
#> 2 0.14551661  8000 1164.133
```
