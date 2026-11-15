# Seasonal ARIMA forecast with a hold-out check

Fits SARIMA\\(p,d,q)(P,D,Q)\_s\\ (by the method of
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md))
to all but the last `h` values of a count series, forecasts those `h`,
and scores the forecasts against them: the mean absolute percentage
error (each error over the larger of the count and 1) and the root mean
square error. The model must difference the series once at least
(`order[2]` or `seasonal[2]` above 0), so it has no mean.

## Usage

``` r
sarfc(y, h = 12, order = c(1, 1, 1), seasonal = c(0, 1, 1, 12))
```

## Arguments

- y:

  The counts, in time order, 36 at least.

- h:

  The values held out.

- order:

  \\(p, d, q)\\.

- seasonal:

  \\(P, D, Q, s)\\, \\s\\ the period (12 for months).

## Value

A list of class `sarfc`: `fc`, `actual`, `coef`, `aic`, `bic` (with the
training length), `mape` (percent), `rmse` and `ntrain`.

## References

Box, G. E. P. and Jenkins, G. M. (1970). *Time Series Analysis:
Forecasting and Control*. Holden-Day.

Gardner, G., Harvey, A. C. and Phillips, G. D. A. (1980). Algorithm AS
154: an algorithm for exact maximum likelihood estimation of
autoregressive-moving average models by means of Kalman filtering.
*Journal of the Royal Statistical Society C*, 29(3), 311-322.

## See also

Other counts over time:
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
[`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md),
[`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md),
[`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md),
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md)

## Examples

``` r
set.seed(1)
y <- rpois(60, 50 + 10 * sin(2 * pi * (1:60) / 12))
sarfc(y, h = 6)[c("mape", "rmse")]
#> $mape
#> [1] 16.91721
#> 
#> $rmse
#> [1] 9.036581
#> 
```
