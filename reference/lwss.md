# Lowess smoothing

Robust locally weighted regression (Cleveland 1979) with the conventions
of the classic implementation, so that the fit equals that of
[`stats::lowess()`](https://rdrr.io/r/stats/lowess.html); the fitted
values come back in the order of the input, not sorted.

## Usage

``` r
lwss(x, y, f = 2/3, iter = 3, delta = NULL)
```

## Arguments

- x, y:

  The points: finite numbers of equal length.

- f:

  The share of points in each neighbourhood: positive.

- iter:

  The number of robustness steps: a whole number from 0.

- delta:

  Points within `delta` of the last fitted one are interpolated: `NULL`
  for 1 per cent of the range of `x`, or a number from 0.

## Value

A numeric vector of fitted values, one per point in input order.

## References

Cleveland, W. S. (1979). Robust locally weighted regression and
smoothing scatterplots. *Journal of the American Statistical
Association*, 74(368), 829-836.
[doi:10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038)

## Examples

``` r
lwss(c(3, 1, 2, 5, 4), c(2, 1, 3, 6, 4))
#> [1] 2 1 3 6 4
```
