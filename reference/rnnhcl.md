# Risk-adjusted nearest-neighbour hierarchical clustering

Nearest-neighbour clustering (as
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md),
first order) with the threshold distance adjusted for a baseline
population at risk, so that incidents dense only because people are
dense there are not called a cluster. In
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)
the threshold is the mean nearest-neighbour distance under complete
spatial randomness plus `z` standard errors,
`0.5 sqrt(A / n) + z sqrt((4 - pi) A / (4 pi)) / n`, which is
`f / sqrt(lambda)` with `lambda = n / A` and
`f = 0.5 + z sqrt((4 - pi) / (4 pi n))`. Here the uniform intensity is
replaced by the intensity the `n` incidents would have if they followed
the baseline: `lambda_i = n b(p_i) / B`, where `b` is the normal-kernel
density (bandwidth `h`) of the baseline units with sizes `bw` and `B` is
their total, and two incidents `i` and `j` are linked when
`|p_i - p_j| < f / sqrt((lambda_i + lambda_j) / 2)`. Taking the mean of
the two intensities is this function's reading of "the threshold at the
pair's baseline density" (the published description does not say how the
two ends combine); a uniform baseline gives back
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)'s
threshold. Grouping is
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)'s
first order; higher orders are not formed.

## Usage

``` r
rnnhcl(
  x,
  y = NULL,
  base,
  bw = NULL,
  h,
  p = 0.5,
  minpts = 10,
  k = 1,
  c = 1,
  lonlat = FALSE
)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  names them, or `lon` and `lat` with `lonlat = TRUE`) and `y` `NULL`;
  or `x` an sf or sfc object of POINT geometries, whose coordinate
  reference system must agree with `lonlat`. An unnamed two-column
  matrix is read as `x` then `y`; other columns are never taken by
  position. Vectors given separately (`x`, `y`, and `w` or `t` where a
  function takes them) are matched by name when every one is named and
  no name repeats, as elsewhere in the package; otherwise by position.
  No missing or infinite values.

- base:

  The baseline units: coordinates as for `x` (zone centroids, or
  points), in the same system as `x` (degrees when `lonlat = TRUE`).

- bw:

  `NULL` for size 1 each, or the baseline sizes (populations), numbers
  from 0 with a positive sum.

- h:

  The bandwidth of the baseline's normal kernel, positive, in the unit
  of the coordinates (metres when `lonlat = TRUE`).

- p:

  The one-tailed probability that a pair is as close by chance: one
  number strictly between 0 and 1.

- minpts:

  The fewest points a cluster keeps: a whole number from 2.

- k:

  The divisor offset: 0, 1 or 2; the total weight must exceed it.

- c:

  The multiplier of each cluster's standard deviational ellipse (see
  [`sdell()`](https://rootcoder007.github.io/crimlab/reference/sdell.md)).

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list: `cl`, each point's cluster (a row of `clusters`) or `NA`;
`clusters`, a data frame of `points`, `mx`, `my`, `wx`, `wy`, `a`, `b`
and `angle` as in
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md);
`thr`, each point's own threshold `f / sqrt(lambda_i)`.

## References

Johnson, S. C. (1967). Hierarchical clustering schemes. *Psychometrika*,
32(3), 241-254.
[doi:10.1007/BF02289588](https://doi.org/10.1007/BF02289588)

Clark, P. J. and Evans, F. C. (1954). Distance to nearest neighbor as a
measure of spatial relationships in populations. *Ecology*, 35(4),
445-453. [doi:10.2307/1931034](https://doi.org/10.2307/1931034)

## Examples

``` r
set.seed(1)
town <- cbind(x = rnorm(400, 0, 1), y = rnorm(400, 0, 1))
inc <- rbind(town[1:60, ], cbind(x = rnorm(15, 3, 0.1),
                                 y = rnorm(15, 3, 0.1)))
r <- rnnhcl(inc, base = town, h = 0.5, minpts = 5)
r$clusters
#>   points       mx      my      wx       wy         a         b    angle
#> 1     15 2.973022 3.00705 2.97886 3.000247 0.1228748 0.1070088 1.886728
```
