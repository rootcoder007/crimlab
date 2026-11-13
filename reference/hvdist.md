# Great-circle distance

The distance along the sphere between pairs of points given in longitude
and latitude. The central angle is computed as the `atan2` of the cross
and dot products of the two unit vectors, written so that nothing
cancels: accurate to a few units in the last place from coincident
points to antipodes (the spherical case of Vincenty 1975). The Earth is
taken as a sphere of mean radius `r`; distances on the ellipsoid differ
by up to about 0.5 per cent.

## Usage

``` r
hvdist(from, to, r = 6371008.8)
```

## Arguments

- from, to:

  The two ends of each pair, in degrees: a matrix or data frame with
  columns `lon` and `lat` (an unnamed two-column matrix is read as
  longitude then latitude), or one point as `c(lon, lat)`. Both have the
  same number of points, or one of them has one point, which is paired
  with every point of the other. Longitudes from -180 to 180, latitudes
  from -90 to 90; no missing values.

- r:

  The radius of the sphere, positive; the distances come back in its
  unit. The default is the mean Earth radius in metres.

## Value

A numeric vector of distances, one per pair, in the order of the input.

## References

Vincenty, T. (1975). Direct and inverse solutions of geodesics on the
ellipsoid with application of nested equations. *Survey Review*,
23(176), 88-93.
[doi:10.1179/sre.1975.23.176.88](https://doi.org/10.1179/sre.1975.23.176.88)

## See also

[`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
to project longitude and latitude to metres.

## Examples

``` r
# Toronto to Montreal, metres
hvdist(c(-79.38, 43.65), c(-73.57, 45.50))
#> [1] 503922
pts <- data.frame(lon = c(0, 90, 180), lat = c(0, 0, 0))
hvdist(pts, c(0, 0))
#> [1]        0 10007557 20015114
```
