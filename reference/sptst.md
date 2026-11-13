# Andresen's spatial point pattern test

Compares where two point patterns fall over the same areal units
(Andresen 2009): for each unit, the difference between its share of the
test pattern `b` and its share of the base pattern `a`, with a bootstrap
p-value from multinomial resamples of each pattern; the S index is the
share of units whose difference is not significant (1: the patterns are
alike everywhere).

## Usage

``` r
sptst(a, b, nboot = 200, seed = NULL, alpha = 0.05, scale = 100)
```

## Arguments

- a, b:

  The counts of the base and the test pattern by unit: whole numbers
  from 0; `NA` in either leaves the unit out (counted, never taken as
  0).

- nboot:

  The number of bootstrap resamples: a whole number from 1.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

- alpha:

  The level of each unit's test, strictly between 0 and 1.

- scale:

  The multiplier of the differences (100: percentage points).

## Value

A list: `units`, a data frame with one row per unit of `diff`, `p` and
`flag` (the sign of a significant difference, 0 otherwise; `NaN` without
a p-value); `glob`, a named vector: units `left` out, `empty` in both,
`elig`ible and `flagged`, the S index `s`, its `z` against `1 - alpha`,
the exact `p` of so many flags and `nboot`.

## References

Andresen, M. A. (2009). Testing for similarity in area-based spatial
patterns: a nonparametric Monte Carlo approach. *Applied Geography*,
29(3), 333-345.
[doi:10.1016/j.apgeog.2008.12.004](https://doi.org/10.1016/j.apgeog.2008.12.004)

## Examples

``` r
set.seed(1)
a <- rpois(30, 10)
b <- rpois(30, 10)
sptst(a, b, nboot = 199, seed = 1)$glob
#>        left       empty        elig     flagged           s           z 
#>   0.0000000   0.0000000  30.0000000   2.0000000   0.9333333  -0.4188539 
#>           p       nboot 
#>   0.4464579 199.0000000 
```
