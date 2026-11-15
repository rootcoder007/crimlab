# Age-period-cohort model by the intrinsic estimator

Asks how much of the change in a rate (arrests, imprisonment,
victimisation per head) over the years belongs to age, to the period and
to the generation (cohort) people were born into. From a table of counts
by age group (rows) and period (columns), with equal widths, and the
population at risk in each cell, the Poisson log-linear model \$\$\log
E(y\_{ij}) = \log(expo\_{ij}) + \mu + \alpha_i + \beta_j + \gamma\_{j -
i + A}\$\$ has effects of age, period and cohort, each summing to 0.
Because cohort = period - age, the three linear trends cannot be told
apart: the design has one null direction, and the data fit equally well
along it. The intrinsic estimator (Yang, Fu and Land 2004) picks the
solution with no component along that direction, the one a
principal-components (minimum-norm) fit gives; its covariance is the
model-based one projected the same way. The model is fitted by the
iteratively reweighted least squares of
[`glm()`](https://rdrr.io/r/stats/glm.html) (as
[`cntreg()`](https://rootcoder007.github.io/crimlab/reference/cntreg.md)
fits).

## Usage

``` r
apcie(y, expo, level = 0.95, maxit = 100)
```

## Arguments

- y:

  The counts: a numeric matrix with one row per age group (youngest
  first) and one column per period (earliest first), at least three of
  each, numbers from 0, not all 0. Row and column names label the
  effects.

- expo:

  The population at risk of each cell: a matrix of the same dimensions,
  positive.

- level:

  The confidence level of the Wald intervals: one number strictly
  between 0 and 1.

- maxit:

  The steps of the fit: one whole number from 1.

## Value

An object of class `apcie`, a list: `coef`, one row for the intercept
and for every age group, period and cohort (`term`, `level`, `est`,
`se`, `lo`, `hi`; cohorts numbered from 1, the oldest, the oldest age
group in the first period); `fit`, a named vector (`deviance`, `df`,
`loglik`, `iter`, `conv`); and `dims`.

## Details

The constraint is a choice, like any other that identifies the model:
the age, period and cohort trends it reports are one of many that fit
the same counts, and differ from another constraint's by a linear trend
in each. The curvatures (deviations from the trends) are identified
whatever the constraint.

## References

Yang, Y., Fu, W. J. and Land, K. C. (2004). A methodological comparison
of age-period-cohort models: the intrinsic estimator and conventional
generalized linear models. *Sociological Methodology*, 34(1), 75-110.
[doi:10.1111/j.0081-1750.2004.00148.x](https://doi.org/10.1111/j.0081-1750.2004.00148.x)

## Examples

``` r
# arrests by five-year age group (rows) and five-year period (columns)
set.seed(11)
ag <- c(-1, 0.5, 0.2, -0.3, -0.6)
pe <- c(-0.2, 0, 0.1, 0.1)
co <- seq(-0.4, 0.4, length.out = 8)
ck <- outer(1:5, 1:4, function(i, j) j - i + 5)
eta <- outer(ag, pe, "+") + matrix(co[ck], 5)
pop <- matrix(1e5, 5, 4)
y <- matrix(rpois(20, pop * exp(-6 + eta)), 5,
            dimnames = list(c("15-19", "20-24", "25-29", "30-34", "35-39"),
                            c(2000, 2005, 2010, 2015)))
apcie(y, pop)
#> apcie: 5 ages, 4 periods, 8 cohorts; deviance 7.1182 on 6 df
#>         term level           est         se           lo           hi
#>  (Intercept)       -6.2524461958 0.02235638 -6.296263896 -6.208628495
#>          age 15-19 -0.5967910284 0.04565097 -0.686265286 -0.507316770
#>          age 20-24  0.8181347740 0.02738804  0.764455203  0.871814345
#>          age 25-29  0.4571013806 0.03144387  0.395472536  0.518730225
#>          age 30-34 -0.1427156034 0.03780542 -0.216812864 -0.068618343
#>          age 35-39 -0.5357295228 0.04059046 -0.615285353 -0.456173692
#>       period  2000 -0.3837793842 0.03014061 -0.442853889 -0.324704880
#>       period  2005 -0.0003931638 0.02913733 -0.057501287  0.056714960
#>       period  2010  0.1501259914 0.02581925  0.099521189  0.200730794
#>       period  2015  0.2340465566 0.02873443  0.177728103  0.290365010
#>       cohort     1  0.1133769341 0.07568888 -0.034970536  0.261724405
#>       cohort     2 -0.1340621684 0.06597129 -0.263363521 -0.004760816
#>       cohort     3 -0.0913148531 0.05301973 -0.195231623  0.012601917
#>       cohort     4 -0.0613495318 0.04112454 -0.141952150  0.019253087
#>       cohort     5 -0.0350107914 0.03827002 -0.110018660  0.039997077
#>       cohort     6  0.0636191215 0.03380790 -0.002643152  0.129881395
#>       cohort     7  0.0937161963 0.03550707  0.024123611  0.163308782
#>       cohort     8  0.0510250929 0.08760061 -0.120668942  0.222719128
```
