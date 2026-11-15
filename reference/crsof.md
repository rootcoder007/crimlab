# The coordinate reference system of a spatial data object

Gives the system as its EPSG code and a label in the style of WKT2,
never as a PROJ string.

## Usage

``` r
crsof(x)
```

## Arguments

- x:

  A
  [`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md)
  object.

## Value

A list: `epsg` and `wkt`.

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md),
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
crsof(spobj(data.frame(lon = 1, lat = 1)))
#> $epsg
#> [1] 4326
#> 
#> $wkt
#> [1] "GEOGCRS[\"WGS 84\", EPSG:4326]"
#> 
```
