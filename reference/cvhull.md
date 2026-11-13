# Convex hull of points

The smallest convex polygon holding every point, by Andrew's monotone
chain with an exact orientation predicate (Shewchuk 1997), so
near-collinear points never make a wrong turn. Points on an edge are not
vertices; of coincident points the first is reported.

## Usage

``` r
cvhull(x, y = NULL, lonlat = FALSE)
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

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list: `idx`, the row numbers of the vertices, anticlockwise from the
point with the lowest `x` (then lowest `y`); `area` and `perim`, the
hull's area and perimeter (in the square and the unit of the
coordinates, or square metres and metres when `lonlat = TRUE`). One
distinct point gives one vertex and zeros; collinear points give their
two ends, area 0 and twice the segment's length. Each coordinate must be
0 or of magnitude from 1e-140 to 1e150, where the predicate is exact.

## References

Andrew, A. M. (1979). Another efficient algorithm for convex hulls in
two dimensions. *Information Processing Letters*, 9(5), 216-219.
[doi:10.1016/0020-0190(79)90072-3](https://doi.org/10.1016/0020-0190%2879%2990072-3)

Shewchuk, J. R. (1997). Adaptive precision floating-point arithmetic and
fast robust geometric predicates. *Discrete and Computational Geometry*,
18(3), 305-363.
[doi:10.1007/PL00009321](https://doi.org/10.1007/PL00009321)

## Examples

``` r
cvhull(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5))
#> $idx
#> [1] 1 2 3 4
#> 
#> $area
#> [1] 1
#> 
#> $perim
#> [1] 4
#> 
```
