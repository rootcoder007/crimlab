# Local Moran (LISA)

Anselin's (1995) local indicator of spatial association for each unit:
with `z = x - mean(x)` and `m2 = mean(z^2)`,
`I_i = z_i / m2 sum_j w_ij z_j`, its expectation and variance under
total randomisation, the two-sided normal p-value, and a conditional
permutation p-value: unit `i` keeps its value and its neighbours' values
are drawn without replacement from the others, `nsim` times from `seed`.
The quadrant tells high-high and low-low clusters from high-low and
low-high outliers. The p-values are per unit and not adjusted for the
many tests made (use
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)); with
row-standardised weights the sum of the `I_i` is `S0` times global
Moran's I. That identity is proved in Lean: see `proof("lmorn")`.

## Usage

``` r
lmorn(x, w, nsim = 0, seed = NULL)
```

## Arguments

- x:

  The values: finite numbers, one per unit.

- w:

  The weights:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, or an spdep `listw` object, for as many units as `x`.

- nsim:

  The number of permutations: a whole number from 0.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

A data frame with one row per unit: `ii` (I_i), `e`, `v`, `z`, `p`,
`pperm` (`(1 + min(#{I* >= I}, #{I* <= I})) / (nsim + 1)`, `NaN` without
permutations) and `quad` (1 high-high, 2 low-low, 3 low-high, 4
high-low, 0 neither).

## References

Anselin, L. (1995). Local indicators of spatial association: LISA.
*Geographical Analysis*, 27(2), 93-115.
[doi:10.1111/j.1538-4632.1995.tb00338.x](https://doi.org/10.1111/j.1538-4632.1995.tb00338.x)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- ifelse(p[, 1] > 0.7 & p[, 2] > 0.7, 5, 1) + rnorm(60, sd = 0.3)
l <- lmorn(v, spwt(p, k = 5), nsim = 199, seed = 1)
head(l)
#>            ii           e         v          z         p pperm quad
#> 1 -0.08136911 -0.01694915 0.1668323 -0.1577177 0.8746793 0.290    4
#> 2  0.12229676 -0.01694915 0.1668323  0.3409120 0.7331698 0.215    2
#> 3  0.03723394 -0.01694915 0.1668323  0.1326550 0.8944662 0.480    2
#> 4  0.14258312 -0.01694915 0.1668323  0.3905786 0.6961088 0.095    2
#> 5  0.18202573 -0.01694915 0.1668323  0.4871449 0.6261557 0.255    2
#> 6  0.09264972 -0.01694915 0.1668323  0.2683280 0.7884469 0.390    2
```
