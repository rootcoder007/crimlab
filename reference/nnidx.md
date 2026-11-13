# Nearest-neighbour index

The mean distance from each point to its k-th nearest neighbour, against
its expectation under complete spatial randomness (CSR) in the rectangle
`win`: the nearest-neighbour index of Clark and Evans (1954) and its
k-order extension (Thompson 1956), with Donnelly's (1978) edge-corrected
expectation for the first order and, when `nsim > 0`, a Monte Carlo test
against `nsim` patterns of as many uniform points in the rectangle. The
Monte Carlo p-values have exact size; the z statistics ignore that
nearest-neighbour distances are dependent and, uncorrected, the edge.
The index tests only the k-th neighbour distances, not complete spatial
randomness as a whole.

## Usage

``` r
nnidx(
  x,
  y = NULL,
  win,
  kmax = 1,
  nsim = 0,
  seed = NULL,
  dist = FALSE,
  metric = "euclid"
)
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

- kmax:

  The highest neighbour order: one whole number from 1 to the number of
  points less 1.

- nsim:

  The number of simulated patterns for the Monte Carlo test: one whole
  number from 0 (no test).

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

- dist:

  `TRUE` to return every point's neighbour distances.

- metric:

  The distance: `"euclid"` or `"manhattan"` (`|dx| + |dy|`, a street
  grid aligned with the axes), as in
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)
  and
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md).
  Under Manhattan distance the expectations use the area of the
  Manhattan ball, `2 r^2`: the mean k-th distance is `sqrt(pi / 2)`
  times the Euclidean expectation and `E R_k^2 = k / (2 lambda)`; no
  edge-corrected (Donnelly) expectation is published for it, so `ce`
  holds `NaN` except `z`, about the uncorrected expectation with
  `Var R_1 = (4 - pi) / (8 lambda)`. The Monte Carlo test is unchanged.

## Value

A list. `order`, a data frame with one row per order `k`: `mean` (the
mean k-th neighbour distance), `exp` (its expectation under CSR without
edge effects, `c_k / sqrt(lambda)`), `ratio`, `pclu` and `preg` (the
Monte Carlo p-values against clustering and against regularity, `NaN`
without simulations) and `sd` (the standard deviation of one k-th
distance under CSR). `ce`, a named vector for the first order:
`donnelly` (the edge-corrected expectation), `ratio` (the index against
it), `z` and `zc` (the Clark-Evans z statistics about the uncorrected
and the corrected expectation). `dist`, `NULL` or a matrix with one row
per point (input order) and one column per order. Distances are in the
unit of the coordinates.

## References

Clark, P. J. and Evans, F. C. (1954). Distance to nearest neighbor as a
measure of spatial relationships in populations. *Ecology*, 35(4),
445-453. [doi:10.2307/1931034](https://doi.org/10.2307/1931034)

Thompson, H. R. (1956). Distribution of distance to nth neighbour in a
population of randomly distributed individuals. *Ecology*, 37(2),
391-394. [doi:10.2307/1933159](https://doi.org/10.2307/1933159)

Donnelly, K. P. (1978). Simulations to determine the variance and edge
effect of total nearest-neighbour distance. In I. Hodder (ed.),
*Simulation Methods in Archaeology*, 91-95. Cambridge University Press.

## Examples

``` r
set.seed(1)
x <- runif(200)
y <- runif(200)
nnidx(x, y, win = c(0, 1, 0, 1), kmax = 3, nsim = 99)
#> $order
#>   k       mean        exp    ratio pclu preg         sd
#> 1 1 0.03753432 0.03535534 1.061631 0.84 0.17 0.01848106
#> 2 2 0.05462585 0.05303301 1.030035 0.44 0.57 0.01925094
#> 3 3 0.06778844 0.06629126 1.022585 0.19 0.82 0.01949659
#> 
#> $ce
#>  donnelly     ratio         z        zc 
#> 0.0364416 1.0299854 1.6674068 0.8361710 
#> 
#> $dist
#> NULL
#> 
```
