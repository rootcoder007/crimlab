# K-means partition of points

Partitions the points into `k` groups by Lloyd's iterations: each point
joins its nearest centre (the lowest-numbered on ties), then each centre
moves to its points' mean (`upd = "mean"`, K-means, minimising the
within-group sum of squares) or to their centre of minimum distance
(`upd = "weber"`, a spatial K-medians minimising the sum of distances),
until no point changes group or after `maxit` moves. An emptied group
keeps its centre. The result depends on the start: k-means++ seeding
(Arthur and Vassilvitskii 2007), the centres of the busiest grid cells
at least a given distance apart, or centres given.

## Usage

``` r
kmpart(
  x,
  y = NULL,
  k,
  init = "pp",
  cell = NULL,
  upd = "mean",
  maxit = 100,
  seed = NULL,
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

- k:

  The number of groups: one whole number from 1 to the number of points.

- init:

  The start: `"pp"` (k-means++, drawn from `seed`), `"cells"` (square
  cells of side `cell[1]` from the lower-left corner of the points'
  bounding box, taken by decreasing count, a cell's centre becoming a
  seed when it is at least `cell[2]` from every earlier one), or the `k`
  starting centres as a matrix or data frame of coordinates (degrees
  when `lonlat = TRUE`).

- cell:

  For `init = "cells"`: `c(side, separation)`, the side positive and the
  separation from 0, in the unit of the coordinates.

- upd:

  The centre update: `"mean"` or `"weber"`.

- maxit:

  The largest number of iterations: one whole number from 0.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list. `cl`, the group of each point (1 to `k`, input order). `cen`, a
data frame with one row per group: `x`, `y` (its centre, in degrees when
`lonlat = TRUE`), `n` (points), `ss` (sum of squared distances to the
centre) and `sd` (sum of distances). `fit`, a named vector: `ss` and
`sd` (totals), `iter` (moves made) and `conv` (1 if no point changed
group at the end).

## References

Lloyd, S. (1982). Least squares quantization in PCM. *IEEE Transactions
on Information Theory*, 28(2), 129-137.
[doi:10.1109/TIT.1982.1056489](https://doi.org/10.1109/TIT.1982.1056489)

Arthur, D. and Vassilvitskii, S. (2007). k-means++: the advantages of
careful seeding. *Proceedings of the Eighteenth Annual ACM-SIAM
Symposium on Discrete Algorithms*, 1027-1035.

## Examples

``` r
set.seed(1)
x <- c(rnorm(30, 0), rnorm(30, 5))
y <- c(rnorm(30, 0), rnorm(30, 5))
kmpart(x, y, k = 2, seed = 1)$cen
#>            x        y  n       ss       sd
#> 1 5.13277458 5.113333 30 40.91834 31.93770
#> 2 0.08245817 0.110278 30 51.51169 33.83318
```
