# The specification of a web map

Describes an interactive map of the points as data, for a mapping
library to draw: the tile provider, the centre (the mean longitude and
latitude), a zoom level and the points in longitude and latitude. Web
Mercator coordinates (EPSG:3857) are taken back to longitude and
latitude by the inverse of the spherical Mercator projection.

## Usage

``` r
lfspec(x)
```

## Arguments

- x:

  A
  [`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md)
  object.

## Value

A list: `provider`, `interactive`, `center` (`lon`, `lat`), `zoom` and
`points` (a matrix of `lon` and `lat`).

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md),
[`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md),
[`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md),
[`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md),
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md),
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md),
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
lfspec(spobj(data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7))))
#> $provider
#> [1] "OpenStreetMap"
#> 
#> $interactive
#> [1] TRUE
#> 
#> $center
#>    lon    lat 
#> -79.35  43.65 
#> 
#> $zoom
#> [1] 8
#> 
#> $points
#>        lon  lat
#> [1,] -79.4 43.6
#> [2,] -79.3 43.7
#> 
```
