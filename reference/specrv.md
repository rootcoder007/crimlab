# Specification curve

The estimate of a treatment effect under every reasonable specification
at once, and whether the specifications as a whole show an effect: the
specification curve of Simonsohn, Simmons and Nelson (2020). A finding
that holds under one specification only is fragile; one that survives
the choices an analyst could defensibly have made is not.

## Usage

``` r
specrv(
  data,
  y,
  d,
  x = NULL,
  opt = NULL,
  se = c("hc1", "iid"),
  level = 0.95,
  null = c("perm", "boot"),
  nsim = 500,
  seed = NULL
)

# S3 method for class 'specrv'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  The names of the outcome columns (numbers): each is a specification of
  the outcome.

- d:

  The name of the treatment column (numbers; it must vary).

- x:

  `NULL` or the names of the covariates in every specification
  (numbers); for [`print()`](https://rdrr.io/r/base/print.html), an
  object of class `specrv`.

- opt:

  `NULL` or the names of the optional covariates (numbers), at most ten.

- se:

  `"hc1"` or `"iid"`.

- level:

  The level of the intervals and of significance: one number strictly
  between 0 and 1.

- null:

  `"perm"` or `"boot"`.

- nsim:

  The data sets under the null: a whole number from 1 to 100000.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `specrv`, a list: `spec`, one row per specification
(`outcome`, `ctrl`, the optional covariates in it or `"none"`, `est`,
`se`, `stat`, `p`, `lo`, `hi`, `sig`); `test`, one row per summary
(`stat`: `median`, `nsig`, `z`; `obs`; `p`); `sign`, the dominant
direction; `dist`, the summaries under the null (a matrix, `NA` for a
draw left out); `nok`; and `null`, `nsim`, `seed`.

## Details

The specifications are every outcome in `y` crossed with every subset of
the optional covariates `opt` (`2^length(opt)` subsets, from none to
all), the covariates `x` always in. Each is the least squares fit of the
outcome on an intercept, the covariates and `d`, giving the coefficient
of `d` with a classical or HC1 standard error and Student's t. The curve
is summarised three ways (their sec. 3): the median estimate; the number
of specifications significant at `1 - level` in the dominant direction
(the sign of the median); and the sum of the specifications' t
statistics over the square root of their number. Each summary is
compared with its values on `nsim` data sets where the null holds: with
`null = "perm"` the treatment is permuted across rows (for a randomised
treatment); with `null = "boot"` each specification's estimate is
removed from its outcome (`y - b d`) and rows are drawn with
replacement. The p-value of a summary is `(1 + k) / (1 + nok)`, `k` of
the `nok` usable data sets at least as extreme in the dominant
direction; a bootstrap draw whose fit is singular is left out.

## References

Simonsohn, U., Simmons, J. P. and Nelson, L. D. (2020). Specification
curve analysis. *Nature Human Behaviour*, 4(11), 1208-1214.
[doi:10.1038/s41562-020-0912-z](https://doi.org/10.1038/s41562-020-0912-z)

Steegen, S., Tuerlinckx, F., Gelman, A. and Vanpaemel, W. (2016).
Increasing transparency through a multiverse analysis. *Perspectives on
Psychological Science*, 11(5), 702-712.
[doi:10.1177/1745691616658637](https://doi.org/10.1177/1745691616658637)

## Examples

``` r
set.seed(7)
n <- 200
d <- data.frame(t = rbinom(n, 1, 0.5), age = runif(n, 18, 60),
                prior = rpois(n, 2), male = rbinom(n, 1, 0.5))
d$arrest <- 0.3 * d$t + 0.1 * d$prior + rnorm(n)
d$charge <- 0.2 * d$t + 0.1 * d$prior + rnorm(n)
s <- specrv(d, c("arrest", "charge"), "t", opt = c("age", "prior",
            "male"), nsim = 200, seed = 1)
s
#> specrv: 16 specifications, 0 significant; perm null, 200 data sets
#>    stat    obs       p
#>  median 0.1427 0.06468
#>    nsig 0.0000 1.00000
#>       z 3.9272 0.07463
head(s$spec)
#>   outcome        ctrl       est        se     stat          p          lo
#> 1  arrest        none 0.2145197 0.1434022 1.495931 0.13626394 -0.06827189
#> 2  arrest         age 0.2124390 0.1436782 1.478575 0.14085117 -0.07090585
#> 3  arrest       prior 0.2389817 0.1435041 1.665330 0.09743638 -0.04401974
#> 4  arrest age + prior 0.2373957 0.1437597 1.651337 0.10027101 -0.04611867
#> 5  arrest        male 0.2095342 0.1439181 1.455926 0.14700476 -0.07428369
#> 6  arrest  age + male 0.2077787 0.1442450 1.440457 0.15133357 -0.07669274
#>          hi   sig
#> 1 0.4973113 FALSE
#> 2 0.4957839 FALSE
#> 3 0.5219832 FALSE
#> 4 0.5209100 FALSE
#> 5 0.4933520 FALSE
#> 6 0.4922501 FALSE
```
