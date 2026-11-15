# A spatial data object

Holds points with their attributes and their coordinate reference
system: the coordinates come from two columns of a data frame, the other
columns are kept as attributes. The system is an EPSG code, never a PROJ
string: 4326 (WGS 84 longitude and latitude, in degrees), 3857 (Web
Mercator, in metres) or another code, taken to be in metres.
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md)
returns the coordinates,
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md)
the system,
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md)
the specification of a web map and
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md)
a sample of the points.

## Usage

``` r
spobj(data, coords = c("lon", "lat"), crs = 4326)
```

## Arguments

- data:

  A data frame holding the coordinate columns and any attributes.

- coords:

  The names of the two coordinate columns: x (or longitude) first.

- crs:

  The EPSG code: one whole number.

## Value

A list of class `spobj`: `coords` (a two-column matrix), `data` (the
other columns), `cnames`, `crs`, `units` and `n`.

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md),
[`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md),
[`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md),
[`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md),
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md),
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md),
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
d <- data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7), v = c(1, 2))
s <- spobj(d, c("lon", "lat"))
s
#> spobj: 2 points, EPSG:4326 (degrees)
format(s)
#> [1] "spobj: 2 points, EPSG:4326 (degrees)"
summary(s)
#>   coord   min   max
#> 1   lon -79.4 -79.3
#> 2   lat  43.6  43.7
```
