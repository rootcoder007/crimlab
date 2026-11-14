# Mean centre and dispersion of points

The (weighted) mean centre of a point pattern, the standard deviations
of the two coordinates, the standard distance deviation (the root mean
squared distance from the mean centre), and the geometric and harmonic
means of the coordinates (Bachi 1959; Neft 1966). The mean centre is the
point that minimises the sum of squared distances to the points. The
dispersions divide by the total weight less `k`: `k = 1` makes the
squared values unbiased for the variances; `k = 2` is the convention of
some crime-mapping texts; `k = 0` divides by the total weight. That the
mean centre minimises the weighted sum of squared distances is proved in
Lean: see `proof("cgstat")`.

## Usage

``` r
cgstat(x, y = NULL, w = NULL, k = 1, lonlat = FALSE)
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

- w:

  `NULL` for equal weights, or one non-negative weight per point (counts
  at a location, say), with a positive sum.

- k:

  The divisor offset: 0, 1 or 2; the total weight must exceed it.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A named numeric vector: `w` (the total weight), `mx`, `my` (the mean
centre, in degrees when `lonlat = TRUE`), `sx`, `sy` (standard
deviations of the coordinates), `sdd` (the standard distance deviation),
`gx`, `gy` (geometric means) and `hx`, `hy` (harmonic means). The last
four are `NaN` for a coordinate with a non-positive value under positive
weight, and for `lonlat = TRUE` (they depend on the origin of the
coordinates). Lengths are in the unit of the coordinates, or metres.

## References

Bachi, R. (1959). Statistical analysis of geographical series. *Bulletin
de l'Institut International de Statistique*, 36(2), 229-240.

Neft, D. S. (1966). *Statistical Analysis for Areal Distributions*.
Regional Science Research Institute, Philadelphia.

## Examples

``` r
set.seed(1)
p <- cbind(x = rnorm(50, 100, 10), y = rnorm(50, 200, 5))
cgstat(p)
#>          w         mx         my         sx         sy        sdd         gx 
#>  50.000000 101.004483 200.586632   8.313939   4.844140   9.622228 100.653418 
#>         gy         hx         hy 
#> 200.529543 100.284986 200.472685 
cgstat(p[, 1], p[, 2], w = rep(1:2, 25))
#>          w         mx         my         sx         sy        sdd         gx 
#>  75.000000 100.685776 200.457192   8.553285   4.753294   9.785321 100.307901 
#>         gy         hx         hy 
#> 200.401917  99.908680 200.346969 
cgstat(p, k = 0)  # dispersions over n, not n - 1
#>          w         mx         my         sx         sy        sdd         gx 
#>  50.000000 101.004483 200.586632   8.230380   4.795453   9.525520 100.653418 
#>         gy         hx         hy 
#> 200.529543 100.284986 200.472685 
```
