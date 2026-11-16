# Least squares with regression diagnostics and predictions

The multiple regression that Algahtany, Kumar and Barclay (2022) fit to
monthly counts of each kind of offence in two regions, with the mean
temperature, the mean relative humidity and the days of haze in the
month as predictors, together with the checks they report before reading
it: standardised coefficients, the tolerance and variance inflation of
each predictor and the correlations among them, the Durbin-Watson
statistic of the residuals, each row's leverage, studentised residuals
and Cook's distance, and the expected count, with its standard error and
interval, for a month of given weather.

## Usage

``` r
olsdx(y, x, at = NULL, int = TRUE, pred = "mean", level = 0.95)

# S3 method for class 'olsdx'
print(x, cook = 1, lev = c(0.2, 0.5), r = 0.7, vif = 10, ...)
```

## Arguments

- y:

  The outcome of each row (a monthly count of offences, say): finite
  numbers, not all equal.

- x:

  The predictors: a numeric matrix or data frame, one row per value of
  `y` (or one vector); column names name the terms. For
  [`print()`](https://rdrr.io/r/base/print.html), an `olsdx` result.

- at:

  `NULL`, or the points at which to predict: a numeric matrix or data
  frame with the columns of `x` in their order, one row per point, or
  one point as a vector.

- int:

  `TRUE` (default) for a constant in the model.

- pred:

  `"mean"` (default) for the interval of the mean at each point, `"new"`
  for that of one new row.

- level:

  The level of the intervals: one number between 0 and 1.

- cook:

  The Cook's distance above which
  [`print()`](https://rdrr.io/r/base/print.html) counts a row.

- lev:

  Two leverages: [`print()`](https://rdrr.io/r/base/print.html) counts
  the rows between them and those above the second.

- r:

  The absolute correlation above which
  [`print()`](https://rdrr.io/r/base/print.html) counts a pair of
  predictors.

- vif:

  The variance inflation factor above which
  [`print()`](https://rdrr.io/r/base/print.html) counts a predictor.

- ...:

  Not used.

## Value

An object of class `olsdx`, a list: `coef`, one row per term (`term`,
`est`, `se`, `t`, `p`, `lo`, `hi`, `beta`, `tol`, `vif`, the last three
`NA` for the constant); `fit`, one row (`n`, `k`, `r` the multiple
correlation, `r2`, `adjr2`, `sigma`, `f`, `df1`, `df2`, `p` the F
test's, `dw`, `rss`); `cases`, one row per value of `y` (`fit`, `resid`,
`zres`, `sres`, `tres`, `hat`, `cook`); `cor`, the correlations of the
columns of `x`; `vcov`, the variance of the coefficients; `pred`, `NULL`
without `at`, else one row per point (`fit`, `se`, `lo`, `hi`); `int`,
`ptype` (`pred`) and `level`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## Details

The fit is ordinary least squares of `y` on a constant (with `int`) and
the columns of `x`; tests and intervals use Student's t on `n - k`
degrees of freedom, `k` the coefficients. `beta` is `b sd(x) / sd(y)`;
`vif` is `1 / (1 - R2_j)`, `R2_j` that of column `j` regressed on a
constant and the other columns, and `tol` its inverse. `dw` is the sum
of the squared differences of successive residuals over their sum of
squares, taken in the order of the rows: give the months in time order.
In `cases`, `zres` is the residual over `s`, `sres` the internally and
`tres` the externally studentised residual, `hat` the leverage `h` and
`cook` Cook's distance `sres^2 h / (k (1 - h))`: the values of
[`stats::rstandard()`](https://rdrr.io/r/stats/influence.measures.html),
[`stats::rstudent()`](https://rdrr.io/r/stats/influence.measures.html),
[`stats::hatvalues()`](https://rdrr.io/r/stats/influence.measures.html)
and
[`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html).
Without the constant, `r2` and `f` are about 0, as
[`stats::summary.lm()`](https://rdrr.io/r/stats/summary.lm.html) has
them.

At each point of `at` (a row of predictor values), `fit` is the fitted
value and `se` the standard error of the mean of rows with those values
(`pred = "mean"`: the estimate the paper gives, with its interval, for a
month of 35 degrees, 74 per cent humidity and 25 days of haze), or of
one new row (`pred = "new"`, which adds the residual variance: a
prediction interval).

The paper's cut-offs are the arguments of
[`print()`](https://rdrr.io/r/base/print.html), which counts the rows
and predictors beyond them: Cook's distance above 1, leverage from 0.2
to 0.5 and above 0.5, correlations between predictors above 0.7 and a
variance inflation above 10 (a tolerance below 0.1). Its checks of
normality are plots of the standardised residuals, `cases$zres`, as
[`stats::qqnorm()`](https://rdrr.io/r/stats/qqnorm.html) draws them.

## References

Algahtany, M., Kumar, L. and Barclay, E. (2022). A tested method for
assessing and predicting weather-crime associations. *Environmental
Science and Pollution Research*, 29, 75013-75030.
[doi:10.1007/s11356-022-20440-6](https://doi.org/10.1007/s11356-022-20440-6)

## Examples

``` r
# 84 months of weather and assaults, simulated
set.seed(4)
w <- data.frame(temp = rnorm(84, 26, 8), hum = runif(84, 10, 60),
                haze = rpois(84, 15))
assault <- 5 + 0.7 * w$temp + 0.15 * w$hum + 0.35 * w$haze +
  rnorm(84, sd = 4)
r <- olsdx(assault, w, at = c(temp = 35, hum = 74, haze = 25))
r
#> olsdx: 84 rows, 4 coefficients, R2 0.6862, adjusted 0.6744
#> F 58.3 on 3 and 80 df, p = 4.43e-20; Durbin-Watson 1.724
#>         term       est         se         t            p         lo        hi
#>  (Intercept) 4.0071452 2.47472213  1.619230 1.093335e-01 -0.9177087 8.9319992
#>         temp 0.6539668 0.06021597 10.860355 2.087182e-17  0.5341332 0.7738004
#>          hum 0.1836688 0.03019435  6.082888 3.838547e-08  0.1235802 0.2437575
#>         haze 0.3983768 0.09932576  4.010811 1.351956e-04  0.2007123 0.5960414
#>       beta       tol      vif
#>         NA        NA       NA
#>  0.6808282 0.9982123 1.001791
#>  0.3815577 0.9970319 1.002977
#>  0.2513618 0.9987934 1.001208
#> rows with Cook's distance above 1: 0; leverage from 0.2 to 0.5: 0, above: 0
#> pairs of predictors correlated above 0.7: 0; predictors with VIF above 10: 0
#> means at `at`:
#>      fit       se       lo      hi
#>  50.4469 1.693266 47.07719 53.8166
r$pred
#>       fit       se       lo      hi
#> 1 50.4469 1.693266 47.07719 53.8166
# one new month rather than the mean of such months
olsdx(assault, w, at = c(temp = 35, hum = 74, haze = 25),
      pred = "new")$pred
#>       fit       se       lo       hi
#> 1 50.4469 4.289069 41.91138 58.98242
# a stricter count of influential rows
print(r, cook = 4 / 84)
#> olsdx: 84 rows, 4 coefficients, R2 0.6862, adjusted 0.6744
#> F 58.3 on 3 and 80 df, p = 4.43e-20; Durbin-Watson 1.724
#>         term       est         se         t            p         lo        hi
#>  (Intercept) 4.0071452 2.47472213  1.619230 1.093335e-01 -0.9177087 8.9319992
#>         temp 0.6539668 0.06021597 10.860355 2.087182e-17  0.5341332 0.7738004
#>          hum 0.1836688 0.03019435  6.082888 3.838547e-08  0.1235802 0.2437575
#>         haze 0.3983768 0.09932576  4.010811 1.351956e-04  0.2007123 0.5960414
#>       beta       tol      vif
#>         NA        NA       NA
#>  0.6808282 0.9982123 1.001791
#>  0.3815577 0.9970319 1.002977
#>  0.2513618 0.9987934 1.001208
#> rows with Cook's distance above 0.04761905: 6; leverage from 0.2 to 0.5: 0, above: 0
#> pairs of predictors correlated above 0.7: 0; predictors with VIF above 10: 0
#> means at `at`:
#>      fit       se       lo      hi
#>  50.4469 1.693266 47.07719 53.8166
```
