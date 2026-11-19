# Simulated compositions from a Dirichlet distribution

Compositions that are closed (each row sums to `total`) and positive,
such as the element shares a portable X-ray fluorescence scan reports,
drawn from the Dirichlet distribution with concentrations `alpha`: each
part is a gamma variate of shape `alpha[j]` and rate 1 divided by the
sum of its row. Normal noise is the wrong model for such data, as it
leaves the simplex. The draws are for testing a compositional analysis
before real measurements exist
([`clrilr()`](https://rootcoder007.github.io/crimlab/reference/clrilr.md),
[`tpbhm()`](https://rootcoder007.github.io/crimlab/reference/tpbhm.md)),
never a substitute for them. The mean share of part `j` is
`alpha[j] / sum(alpha)`; larger concentrations give less spread around
it.

## Usage

``` r
dirsim(n, alpha, total = 1, seed = NULL)
```

## Arguments

- n:

  The number of compositions: a whole number from 1 to 2147483647.

- alpha:

  The concentration of each part: positive finite numbers, at least two;
  their names name the columns.

- total:

  What each composition sums to: one positive number (1 for shares, 1e6
  for parts per million).

- seed:

  `NULL` to draw the seed from R's generator, or one whole number from
  0.

## Value

A numeric matrix of `n` rows and one column per part.

## Details

The gammas are drawn part by part, `n` of each, from R's generator
(Mersenne-Twister with inversion, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) seeds it), so that
`set.seed(seed); g <-` `matrix(rgamma(n * k,`
`rep(alpha, each = n)), n)` gives the same draws and
`g / rowSums(g) * total` the same result. Without `seed` one is drawn
from R's generator.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.
[doi:10.1007/978-94-009-4109-0](https://doi.org/10.1007/978-94-009-4109-0)

## Examples

``` r
# a calcium-rich profile and a natural one, in parts per million
lime <- c(Ca = 85, P = 5, Fe = 2, Sr = 0.5, Pb = 0.2, Zn = 0.2)
soil <- c(Ca = 30, P = 15, Fe = 5, Sr = 1, Pb = 0.5, Zn = 0.5)
x <- rbind(dirsim(5, lime, 1e6, seed = 1), dirsim(5, soil, 1e6, seed = 2))
round(x)
#>           Ca      P     Fe    Sr   Pb    Zn
#>  [1,] 916028  64989   9709  8530  338   406
#>  [2,] 924519  59016  10998  3742 1585   139
#>  [3,] 926699  55681  13866     7  343  3403
#>  [4,] 948999  41626   6514    75    0  2786
#>  [5,] 918903  38403  34643  7958   93     0
#>  [6,] 545725 298947  82563 67504 1473  3788
#>  [7,] 656421 315410  24984   307 2854    24
#>  [8,] 490101 333923 148857 25907   52  1160
#>  [9,] 549398 336330  94906  4794  858 13714
#> [10,] 553314 280049  72463 29418 2497 62259
clrilr(x)[, "Ca"]
#>  [1] 4.582492 4.644030 5.378437 6.646381 6.411516 2.577556 4.556793 3.286704
#>  [9] 2.856785 2.205577
```
