# The general G of Getis and Ord

The global concentration statistic of Getis and Ord (1992) for
non-negative values: `G = sum_ij w_ij x_i x_j / sum_{i != j} x_i x_j`,
its expectation `S0 / (n (n - 1))` and variance under randomisation, the
z statistic, its upper normal tail (high values near high values) and,
with `nsim > 0`, a permutation p-value. Unlike Moran's I and Geary's C
it tells high-value clusters from low-value ones. Its expectation under
randomisation is proved in Lean: see `proof("gorg")`.

## Usage

``` r
gorg(x, w, nsim = 0, seed = NULL)
```

## Arguments

- x:

  The values: finite numbers from 0, at least two positive.

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

A named numeric vector: `g`, `e`, `v`, `z`, `p` and `pperm`
(`(1 + #{G* >= G}) / (nsim + 1)`).

## References

Getis, A. and Ord, J. K. (1992). The analysis of spatial association by
use of distance statistics. *Geographical Analysis*, 24(3), 189-206.
[doi:10.1111/j.1538-4632.1992.tb00261.x](https://doi.org/10.1111/j.1538-4632.1992.tb00261.x)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
gorg(exp(3 * p[, 1]), spwt(p, type = "band", d = 0.2, style = "raw"))
#>            g            e            v            z            p        pperm 
#> 1.663958e-01 1.167347e-01 1.514881e-04 4.034852e+00 2.731842e-05          NaN 
```
