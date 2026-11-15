# The coordinates of a spatial data object

The coordinates of a spatial data object

## Usage

``` r
spxy(x)
```

## Arguments

- x:

  A
  [`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md)
  object.

## Value

The two-column matrix of coordinates, named by their columns.

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md),
[`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md),
[`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md),
[`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md),
[`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md),
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md),
[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md),
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
spxy(spobj(data.frame(lon = 1:3, lat = 1:3)))
#>      lon lat
#> [1,]   1   1
#> [2,]   2   2
#> [3,]   3   3
```
