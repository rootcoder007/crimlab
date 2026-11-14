# Median centre of points

The median of the `x` coordinates and the median of the `y` coordinates
(the mean of the two middle values for an even number of points). It
depends on the orientation of the axes: rotating the points does not
rotate it.
[`webpt()`](https://rootcoder007.github.io/crimlab/reference/webpt.md)
gives the rotation-invariant centre of minimum distance.

## Usage

``` r
cgmedn(x, y = NULL, lonlat = FALSE)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as `sf::st_coordinates()` names them, or `lon` and `lat` with
  `lonlat = TRUE`) and `y` `NULL`; or `x` an sf or sfc object of POINT
  geometries, whose coordinate reference system must agree with
  `lonlat`. An unnamed two-column matrix is read as `x` then `y`; other
  columns are never taken by position. Vectors given separately (`x`,
  `y`, and `w` or `t` where a function takes them) are matched by name
  when every one is named and no name repeats, as elsewhere in the
  package; otherwise by position. No missing or infinite values.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A named numeric vector `c(x, y)`, in degrees when `lonlat = TRUE` (the
medians are taken of the projected coordinates).

## References

Neft, D. S. (1966). *Statistical Analysis for Areal Distributions*.
Regional Science Research Institute, Philadelphia.

## Examples

``` r
cgmedn(c(1, 2, 3, 10), c(5, 1, 2, 2))
#>   x   y 
#> 2.5 2.0 
```
