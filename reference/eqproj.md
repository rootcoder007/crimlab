# Equal-area projection of longitude and latitude, and its inverse

`eqproj()` maps longitude and latitude to planar metres by the Lambert
azimuthal equal-area projection of the sphere about a centre; `eqinv()`
maps them back. Areas are exact, so intensities, nearest-neighbour
expectations and K functions computed on the result are right; lengths
are scaled by `cos(c / 2)` along the radius from the centre and its
inverse across it (`c` the central angle), so within 100 km of the
centre every length is right to 3.1e-5 relative. Every planar function
of the package takes the projected coordinates; those with a `lonlat`
argument call this projection themselves.

## Usage

``` r
eqproj(x, y = NULL, lon0 = NULL, lat0 = NULL, r = 6371008.8)

eqinv(x, y = NULL, lon0, lat0, r = 6371008.8)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors, or `x` a matrix or data frame
  with columns named `lon` and `lat` for `eqproj()` and `x` and `y` for
  `eqinv()` (an unnamed two-column matrix is read in that order) and `y`
  `NULL`. Longitudes from -180 to 180, latitudes from -90 to 90; no
  missing values.

- lon0, lat0:

  The centre of the projection in degrees. For `eqproj()`, `NULL` (both)
  takes the circular mean of the longitudes and the mean of the
  latitudes.

- r:

  The radius of the sphere, positive; the distances come back in its
  unit. The default is the mean Earth radius in metres.

## Value

A data frame with one row per point, in input order: `x` and `y` (east
and north, in the unit of `r`) from `eqproj()`, with the centre used as
the attribute `centre`; `lon` and `lat` from `eqinv()`.

## References

Snyder, J. P. (1987). *Map Projections: A Working Manual*. U.S.
Geological Survey Professional Paper 1395, chapter 24.
[doi:10.3133/pp1395](https://doi.org/10.3133/pp1395)

## Examples

``` r
p <- eqproj(c(-79.40, -79.38, -79.35), c(43.64, 43.66, 43.65))
p
#>            x             y
#> 1 -1877.6518 -1111.6869087
#> 2  -268.1467  1111.9561852
#> 3  2145.5306     0.3446339
a <- attr(p, "centre")
eqinv(p, lon0 = a[1], lat0 = a[2])
#>      lon   lat
#> 1 -79.40 43.64
#> 2 -79.38 43.66
#> 3 -79.35 43.65
```
