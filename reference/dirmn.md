# Directional mean and circular variance

The mean direction of a set of directions, the mean resultant length `R`
(from 0, directions spread evenly, to 1, all the same), the circular
variance `1 - R` and the circular standard deviation `sqrt(-2 log R)`
(Mardia and Jupp 2000). Directions are given as angles, or as points
seen from an origin (the directions of incidents from an offender's
base, say). The mean direction of points depends on the origin chosen.
That `R` lies in `[0, 1]` is proved in Lean: see `proof("dirmn")`.

## Usage

``` r
dirmn(
  x = NULL,
  y = NULL,
  origin = c(0, 0),
  w = NULL,
  angle = NULL,
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

- origin:

  For points: the origin `c(x, y)` (in degrees when `lonlat = TRUE`); no
  point may lie on it.

- w:

  `NULL`, or one weight from 0 per direction, positive sum.

- angle:

  Instead of points: angles in radians, anticlockwise from east (the x
  axis).

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A named numeric vector: `theta` (the mean direction in radians
anticlockwise from east, in `[0, 2 pi)`; `NaN` when `R` is 0), `bearing`
(the same as a compass bearing, degrees clockwise from north), `r`,
`var`, `sd`, `w` (the total weight) and `n`.

## References

Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*. Wiley.
[doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

## Examples

``` r
dirmn(c(1, 2, 0), c(1, 1, 3), origin = c(0, 0))
#>      theta    bearing          r        var         sd          w          n 
#>  0.9315290 36.6273216  0.8948002  0.1051998  0.4714971  3.0000000  3.0000000 
dirmn(angle = c(0.1, 0.2, 6.2))
#>        theta      bearing            r          var           sd            w 
#>  0.072363804 85.853859466  0.993136769  0.006863231  0.117361758  3.000000000 
#>            n 
#>  3.000000000 
```
