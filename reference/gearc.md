# Geary's C

Geary's (1954) contiguity ratio of values on units with spatial weights
`w`: `C = (n - 1) sum w_ij (x_i - x_j)^2 / (2 S0 sum (x_i - xbar)^2)`, 1
under no autocorrelation, below 1 when neighbours are alike. Its
variance is given under normality and under randomisation (Cliff and Ord
1981), with z statistics `(1 - C) / sd` and their upper normal tails,
and, with `nsim > 0`, a permutation p-value against positive
autocorrelation. That C has expectation 1 under randomisation is proved
in Lean: see `proof("gearc")`.

## Usage

``` r
gearc(x, w, nsim = 0, seed = NULL)
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

A named numeric vector: `c`, `vn` and `vr` (its variance under normality
and randomisation), `zn`, `zr`, `pn`, `pr`, `pperm`
(`(1 + #{C* <= C}) / (nsim + 1)`, `NaN` without permutations) and `s0`
(the sum of the weights).

## References

Geary, R. C. (1954). The contiguity ratio and statistical mapping. *The
Incorporated Statistician*, 5(3), 115-145.
[doi:10.2307/2986645](https://doi.org/10.2307/2986645)

Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
Applications*. Pion, London.

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
gearc(p[, 1] + p[, 2] + rnorm(50, sd = 0.1), spwt(p, k = 4))
#>            c           vn           vr           zn           zr           pn 
#> 9.221870e-02 9.592157e-03 9.036261e-03 9.268791e+00 9.549637e+00 9.413378e-21 
#>           pr        pperm           s0 
#> 6.507588e-22          NaN 5.000000e+01 
```
