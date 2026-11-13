# Local Getis-Ord Gi\* on a grid of counts

The local statistic Gi\* of Getis and Ord (1992; Ord and Getis 1995) for
each cell of a grid of counts, with binary weights on the cells within
`rad` cell widths (`rad = 1.5`: the 3 by 3 block, cell included), its
normal p-value and Benjamini-Hochberg q-value, a conditional permutation
p-value, and an exposure-aware version: given the total, the count of a
neighbourhood is binomial with the neighbourhood's share of the exposure
under complete spatial randomness, so cells overhanging the study area
neither inflate the mean nor get flagged. No threshold is applied.

## Usage

``` r
gstr(cnt, ex = NULL, rad = 1.5, nsim = 0, seed = NULL)
```

## Arguments

- cnt:

  Counts per cell: a matrix of whole numbers from 0.

- ex:

  `NULL` for equal exposure, or a matrix of the same dimensions: each
  cell's exposure (its area inside the study window, say), 0 for a cell
  outside (its results are `NaN`).

- rad:

  The neighbourhood radius in cell widths, from 0 to 1000 (1: the cell
  and its four rook neighbours; 1.5: the queen neighbourhood).

- nsim:

  The number of permutations: a whole number from 0.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

A data frame with one row per cell (down the columns of `cnt`): `row`,
`col`, `gz` (the classic Gi\* z), `pz`, `qz`, `pp` (the permutation
p-value, `NaN` without permutations), `ze` (the exposure-aware z), `pe`
(its exact two-sided binomial p-value) and `qe`.

## References

Getis, A. and Ord, J. K. (1992). The analysis of spatial association by
use of distance statistics. *Geographical Analysis*, 24(3), 189-206.
[doi:10.1111/j.1538-4632.1992.tb00261.x](https://doi.org/10.1111/j.1538-4632.1992.tb00261.x)

Ord, J. K. and Getis, A. (1995). Local spatial autocorrelation
statistics: distributional issues and an application. *Geographical
Analysis*, 27(4), 286-306.
[doi:10.1111/j.1538-4632.1995.tb00912.x](https://doi.org/10.1111/j.1538-4632.1995.tb00912.x)

## Examples

``` r
set.seed(1)
m <- matrix(rpois(100, 2), 10)
m[4:5, 4:5] <- 9
g <- gstr(m)
g[g$qe < 0.05, ]
#>    row col       gz           pz           qz  pp       ze           pe
#> 34   4   4 4.374058 1.219581e-05 0.0006097906 NaN 5.368464 2.971884e-06
#> 35   5   4 4.186330 2.835012e-05 0.0009450040 NaN 5.138058 7.201663e-06
#> 44   4   5 4.374058 1.219581e-05 0.0006097906 NaN 5.368464 2.971884e-06
#> 45   5   5 3.998602 6.371769e-05 0.0015929423 NaN 4.907652 1.698365e-05
#>              qe
#> 34 0.0001485942
#> 35 0.0002400554
#> 44 0.0001485942
#> 45 0.0004245913
```
