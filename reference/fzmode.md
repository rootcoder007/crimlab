# Fuzzy mode

For each point, the number of points (itself and coincident ones
included) within distance `r` of it. With `r = 0` it counts coincident
points, the mode of exact locations; with `r > 0` it is the fuzzy mode
of the crime-mapping literature, a fixed-radius scan centred on each
incident. Counts of nearby points overlap by design: each point counts
its own neighbourhood.

## Usage

``` r
fzmode(x, y = NULL, r, lonlat = FALSE)
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

- r:

  The radius: one number from 0, in the unit of the coordinates (metres
  when `lonlat = TRUE`). Distances up to and including `r` count.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

An integer vector, one count per point, in input order.

## References

Eck, J. E., Chainey, S., Cameron, J. G., Leitner, M. and Wilson, R. E.
(2005). *Mapping Crime: Understanding Hot Spots*. National Institute of
Justice, NCJ 209393.

## Examples

``` r
fzmode(c(0, 0, 1, 5), c(0, 0, 0, 5), r = 1)
#> [1] 3 3 3 1
```
