# Global Moran's I with its inference

Moran's (1950) I of values on units with spatial weights `w`: with
`z = x - mean(x)`, `I = n / S0 sum_ij w_ij z_i z_j / sum z_i^2`, its
expectation `-1 / (n - 1)` and its variance under normality and under
randomisation (Cliff and Ord 1981), the z statistics and their upper
normal tails (positive autocorrelation: neighbours alike), and, with
`nsim > 0`, the permutation p-value `(1 + #{I* >= I}) / (nsim + 1)`, the
values shuffled over the units from `seed`. These are the statistics of
[`spdep::moran.test()`](https://r-spatial.github.io/spdep/reference/moran.test.html)
and
[`spdep::moran.mc()`](https://r-spatial.github.io/spdep/reference/moran.mc.html).
Given a matrix or data frame, each column is tested on the same weights:
one row per category or year, a sweep of a variable across a polygon
surface.
[`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md)
gives the local statistics, which add up to `S0 I` with row-standardised
weights.

## Usage

``` r
gmoran(x, w, nsim = 0, seed = NULL)
```

## Arguments

- x:

  The values: finite numbers, one per unit; or a matrix or data frame
  with one column per variable.

- w:

  The weights:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` object, or an `n x n` matrix of weights from 0
  (zero diagonal).

- nsim:

  The number of permutations: a whole number from 0.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

A data frame with one row per variable: `var` (its name), `i`, `e`, `vn`
and `vr`, `zn` and `zr`, `pn` and `pr`, `pperm` (`NaN` without
permutations) and `s0` (the sum of the weights); `NA` but `s0` for a
constant variable.

## References

Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
*Biometrika*, 37(1-2), 17-23.
[doi:10.2307/2332142](https://doi.org/10.2307/2332142)

Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
Applications*. Pion, London.

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
w <- spwt(p, k = 4)
gmoran(p[, 1] + rnorm(50, sd = 0.2), w, nsim = 99, seed = 1)
#>   var         i           e          vn          vr      zn       zr
#> 1   x 0.6251582 -0.02040816 0.007802795 0.007875543 7.30829 7.274458
#>             pn           pr pperm s0
#> 1 1.352813e-13 1.739075e-13  0.01 50
# a sweep over three yearly counts
gmoran(data.frame(y1 = rpois(50, 3), y2 = rpois(50, 3 + 4 * p[, 1]),
                  y3 = rpois(50, 9 * p[, 1])), w)
#>   var             i           e          vn          vr        zn        zr
#> 1  y1 -0.0008085382 -0.02040816 0.007802795 0.007858451 0.2218823 0.2210952
#> 2  y2  0.1023809524 -0.02040816 0.007802795 0.007581203 1.3900640 1.4102328
#> 3  y3  0.3831195986 -0.02040816 0.007802795 0.007584193 4.5682339 4.6336020
#>             pn           pr pperm s0
#> 1 4.122028e-01 4.125092e-01   NaN 50
#> 2 8.225473e-02 7.923548e-02   NaN 50
#> 3 2.459255e-06 1.796788e-06   NaN 50
```
