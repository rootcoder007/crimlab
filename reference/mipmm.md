# Multiple imputation by predictive mean matching

Fills the missing values of numeric columns `m` times by predictive mean
matching, each completed data set from its own chain of chained
equations, and pools an analysis of the completed data sets by Rubin's
rules: the step after
[`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md)
has said how the values are missing. Predictive mean matching imputes
only values that were observed, so an imputed count stays a whole number
and an imputed age stays in its range (van Buuren 2018, sec. 3.4).

## Usage

``` r
mipmm(
  data,
  y,
  x = NULL,
  m = 5,
  k = 5,
  maxit = 5,
  ridge = 1e-05,
  fy = NULL,
  fx = NULL,
  level = 0.95,
  seed = NULL
)

# S3 method for class 'mipmm'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  The names of the numeric columns to impute (missing values as `NA`;
  observed values finite).

- x:

  `NULL` or the names of complete numeric columns that predict them; for
  [`print()`](https://rdrr.io/r/base/print.html), an object of class
  `mipmm`.

- m:

  The imputations: a whole number from 2 to 1000.

- k:

  The donors to draw from: a whole number from 1 to 10000.

- maxit:

  The rounds of each chain: a whole number from 1 to 1000.

- ridge:

  The ridge, as a share of the diagonal: a number from 0.

- fy:

  `NULL` or the outcome of the analysis: a column of `y` or `x`.

- fx:

  `NULL` or its predictors: columns of `y` or `x` other than `fy`.

- level:

  The level of the intervals: one number strictly between 0 and 1.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `mipmm`, a list: `imp`, for each column of `y` a
matrix of its imputed values, one row per missing value (named by its
row of `data`) and one column per imputation; `pool`, a data frame
(`term`, `est`, `se`, `df`, `lo`, `hi`, `p`, `riv`, `fmi`); `nmis`, the
missing values of each column; and `m`, `k`, `maxit`, `seed`.

## Details

Each chain starts the missing values of every column of `y` at draws
from its observed values and, `maxit` times, imputes each column in turn
from all the others (the other columns of `y` as last imputed, and the
complete columns `x`): the least squares fit on its observed rows,
`ridge` times the diagonal of the cross-products added to it; a draw of
the residual scale (`rss` over a chi-square on the residual degrees of
freedom) and of the coefficients about their estimate; and, for each
missing row, a donor drawn with equal chances among the `k` observed
rows whose fitted values (from the estimate) are nearest its own (from
the draw), whose observed value it takes (van Buuren 2018, sec. 4.5;
type 1 matching). The draws come from the package's generator.

The analysis is the least squares fit of `fy` on an intercept and `fx`
(classical variances), or, without `fy`, the mean of each column of `y`
(variance `s^2 / n`). Pooled by Rubin (1987): `est` the mean of the `m`
estimates, `se` the square root of `T = U + (1 + 1 / m) B` (`U` the mean
variance, `B` the variance between imputations), `df` that of Barnard
and Rubin (1999) with the complete-data degrees of freedom of the
analysis, Student's t interval and p-value, `riv` the relative increase
in variance `(1 + 1 / m) B / U` and `fmi` the fraction of missing
information, as the `mice` package pools.

## References

van Buuren, S. (2018). *Flexible Imputation of Missing Data*, 2nd ed.
Chapman and Hall/CRC.
[doi:10.1201/9780429492259](https://doi.org/10.1201/9780429492259)

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
Wiley.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

Barnard, J. and Rubin, D. B. (1999). Small-sample degrees of freedom
with multiple imputation. *Biometrika*, 86(4), 948-955.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

Little, R. J. A. (1988). Missing-data adjustments in large surveys.
*Journal of Business and Economic Statistics*, 6(3), 287-296.
[doi:10.1080/07350015.1988.10509663](https://doi.org/10.1080/07350015.1988.10509663)

## Examples

``` r
set.seed(11)
n <- 200
d <- data.frame(age = round(runif(n, 18, 60)), male = rbinom(n, 1, 0.5))
d$priors <- rpois(n, exp(-1 + 0.04 * d$age))
d$days <- round(20 + 3 * d$priors + 5 * d$male + rnorm(n, sd = 5))
d$priors[sample(n, 40)] <- NA
d$days[sample(n, 30)] <- NA
r <- mipmm(d, c("priors", "days"), c("age", "male"), fy = "days",
           fx = c("priors", "male"), seed = 1)
r
#> mipmm: 70 missing value(s) in 2 column(s), 5 imputations by predictive mean matching
#>         term    est     se    df     lo     hi         p    riv    fmi
#>  (Intercept) 18.629 0.8058 40.95 17.002 20.257 4.065e-25 0.3591 0.2977
#>       priors  3.432 0.2773 27.11  2.863  4.001 1.154e-12 0.5178 0.3849
#>         male  6.242 0.8145 61.78  4.614  7.871 1.550e-10 0.2468 0.2227
head(r$imp$priors)
#>    [,1] [,2] [,3] [,4] [,5]
#> 11    2    1    0    1    0
#> 17    0    1    2    0    2
#> 18    1    0    0    3    2
#> 19    0    0    0    0    2
#> 21    1    1    3    1    1
#> 23    0    1    0    0    0
```
