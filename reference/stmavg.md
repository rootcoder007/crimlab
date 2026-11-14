# Spatial-temporal moving average

For each event, the mean centre of the `span` events centred on it in
time order (ties in input order), the window shortened at the ends of
the series. It traces how the centre of a series of incidents moves over
time. It is a descriptive device of the crime-mapping literature with no
original paper to cite; the mean centre is that of
[`cgstat()`](https://rootcoder007.github.io/crimlab/reference/cgstat.md).

## Usage

``` r
stmavg(x, y = NULL, t, span = 5, lonlat = FALSE)
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

- t:

  The time of each event: finite numbers (for dates,
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) of them, in days
  or seconds).

- span:

  The window: one odd whole number from 1.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A data frame with one row per event, in input order: `mx`, `my` (the
window's mean centre, in degrees when `lonlat = TRUE`) and `n` (events
in the window).

## Examples

``` r
stmavg(x = 1:7, y = c(1, 1, 2, 3, 5, 8, 13), t = 7:1, span = 3)
#>    mx        my n
#> 1 1.5  1.000000 2
#> 2 2.0  1.333333 3
#> 3 3.0  2.000000 3
#> 4 4.0  3.333333 3
#> 5 5.0  5.333333 3
#> 6 6.0  8.666667 3
#> 7 6.5 10.500000 2
```
