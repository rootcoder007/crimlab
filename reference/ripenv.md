# Envelopes of Ripley's L function and a test of complete spatial randomness

`L(r) - r` for the data, its pointwise minimum, maximum and mean over
`nsim` patterns of as many uniform points in the rectangle, and a global
Monte Carlo test on the maximum absolute deviation
`T = max_r |L(r) - r|`, whose p-value `(1 + #{T_sim >= T}) / (nsim + 1)`
has exact size (Ripley 1977; Diggle 2013). The pointwise envelopes are
not a test of the whole curve: each radius alone falls outside them with
probability `2 / (nsim + 1)` under randomness.

## Usage

``` r
ripenv(x, y = NULL, win, r, corr = "iso", nsim = 99, seed = NULL)
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y` (an unnamed two-column
  matrix is read as `x` then `y`). At least two points, all inside
  `win`; no missing values.

- win:

  The study rectangle `c(xmin, xmax, ymin, ymax)`, in the unit of the
  coordinates. The intensity is the number of points over its area.

- r:

  The radii: non-negative numbers in strictly increasing order, in the
  unit of the coordinates.

- corr:

  The edge correction: `"iso"`, `"border"`, `"modif"` or `"none"`.

- nsim:

  The number of simulated patterns: one whole number from 1.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

A list. `env`, a data frame with one row per radius: `r`, `obs`
(`L(r) - r` of the data), `lo`, `hi` (the pointwise envelope) and
`mean`. `test`, a named vector: `t` (the maximum absolute deviation) and
`p` (its Monte Carlo p-value).

## References

Diggle, P. J. (2013). *Statistical Analysis of Spatial and
Spatio-Temporal Point Patterns* (3rd ed.), chapter 6. CRC Press.

Ripley, B. D. (1977). Modelling spatial patterns. *Journal of the Royal
Statistical Society, Series B*, 39(2), 172-212.
[doi:10.1111/j.2517-6161.1977.tb01615.x](https://doi.org/10.1111/j.2517-6161.1977.tb01615.x)

## Examples

``` r
set.seed(1)
ripenv(runif(100), runif(100), win = c(0, 1, 0, 1),
       r = seq(0.01, 0.25, 0.02), nsim = 39)$test
#>          t          p 
#> 0.01033936 0.15000000 
```
