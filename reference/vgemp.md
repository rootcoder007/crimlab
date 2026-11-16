# Empirical semivariogram

The classical estimator of Matheron (1963): half the mean squared
difference of the values over the pairs of points whose distance falls
in each of `nbins` equal bins of `(0, cutoff]` (the first bin closed at
0). Bins without pairs are left out. The semivariogram describes how
dissimilar values become with distance, the first step of kriging: rates
of a crime, a disease or a survey outcome measured at places.

## Usage

``` r
vgemp(x, y = NULL, v, nbins = 15, cutoff = NULL)
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

- v:

  The values: finite numbers, one per point.

- nbins:

  The number of distance bins: one whole number from 1 to 67108864
  (2^26).

- cutoff:

  The largest distance used: `NULL` for a third of the largest distance
  between two points (the convention of gstat), or one positive number.

## Value

A data frame with one row per bin holding pairs: `dist` (the centre of
the bin), `gamma` (the semivariance) and `np` (the number of pairs).

## References

Matheron, G. (1963). Principles of geostatistics. *Economic Geology*,
58(8), 1246-1266.
[doi:10.2113/gsecongeo.58.8.1246](https://doi.org/10.2113/gsecongeo.58.8.1246)

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
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
# simulated values with a trend from west to east
set.seed(1)
p <- cbind(x = runif(60, 0, 10), y = runif(60, 0, 10))
v <- p[, 1] * 0.5 + rnorm(60, 0, 0.3)
vgemp(p, v = v)
#>         dist      gamma np
#> 1  0.3883710 0.03235575 14
#> 2  0.6472850 0.10202814 16
#> 3  0.9061991 0.14740942 25
#> 4  1.1651131 0.17579384 27
#> 5  1.4240271 0.12605216 28
#> 6  1.6829411 0.20833025 43
#> 7  1.9418551 0.27639636 43
#> 8  2.2007691 0.37523216 57
#> 9  2.4596831 0.46234587 48
#> 10 2.7185972 0.54284954 60
#> 11 2.9775112 0.68058945 69
#> 12 3.2364252 0.78158842 69
#> 13 3.4953392 0.74873781 60
#> 14 3.7542532 0.84625342 76
```
