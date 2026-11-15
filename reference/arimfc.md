# ARIMA(1,1,1) forecast of a count series

Fits ARIMA(1,1,1) to a series of counts (monthly incidents, calls for
service, admissions) by exact maximum likelihood, the Kalman filter on
the state-space form with the conditional-sum-of-squares fit as its
start, as [`stats::arima()`](https://rdrr.io/r/stats/arima.html) does by
default, and forecasts `h` periods ahead with their standard errors. AIC
is \\-2 \log L + 2(k + 1)\\; BIC is \\\mathrm{AIC} + (k + 1)(\log n -
2)\\, \\k\\ the ARMA coefficients.

## Usage

``` r
arimfc(y, h = 12)
```

## Arguments

- y:

  The counts, in time order, 24 at least.

- h:

  The periods to forecast.

## Value

A list of class `arimfc`: `fc` and `se` (forecasts and their standard
errors), `coef` (ar1, ma1), `sigma2`, `aic`, `bic` and `n`.

## References

Box, G. E. P. and Jenkins, G. M. (1970). *Time Series Analysis:
Forecasting and Control*. Holden-Day.

Gardner, G., Harvey, A. C. and Phillips, G. D. A. (1980). Algorithm AS
154: an algorithm for exact maximum likelihood estimation of
autoregressive-moving average models by means of Kalman filtering.
*Journal of the Royal Statistical Society C*, 29(3), 311-322.

## See also

Other counts over time:
[`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md),
[`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md),
[`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md),
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
[`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)

## Examples

``` r
# 48 simulated months of counts
set.seed(1)
y <- rpois(48, 50 + 10 * sin(1:48 / 6))
arimfc(y, h = 6)$fc
#> [1] 61.03613 61.69301 61.79991 61.81731 61.82014 61.82060
```
