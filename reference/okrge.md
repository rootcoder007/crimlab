# Ordinary kriging

Predicts the values at new locations `at` by ordinary kriging: the best
linear unbiased predictor under a variogram model with an unknown
constant mean (Cressie 1993, sec. 3.2), with its kriging variance. The
variogram is `vgm`, a list such as
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)
returns; without it, an exponential model is fitted by weighted least
squares to the 15-bin
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md)
semivariogram with the weights `np / h^2` of Cressie (1985), which is
fast but rougher than
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md).

## Usage

``` r
okrge(x, y = NULL, v, at, vgm = NULL)
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

- at:

  The prediction locations: a matrix or data frame of two coordinate
  columns, as `x`.

- vgm:

  `NULL`, or a list with `model` (`"exponential"`, `"spherical"` or
  `"gaussian"`), `nugget`, `psill` (from 0) and `range` (positive).

## Value

A data frame with one row per prediction location: `pred` and `var` (the
kriging variance, at least 0).

## References

Cressie, N. (1985). Fitting variogram models by weighted least squares.
*Mathematical Geology*, 17(5), 563-586.
[doi:10.1007/BF01032109](https://doi.org/10.1007/BF01032109)

Cressie, N. (1993). *Statistics for Spatial Data*, revised edition.
Wiley, New York.
[doi:10.1002/9781119115151](https://doi.org/10.1002/9781119115151)

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md),
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
set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- sin(3 * p[, 1]) + rnorm(60, sd = 0.1)
okrge(p, v = v, at = cbind(0.5, 0.5))
#>      pred        var
#> 1 1.06704 0.01030984
okrge(p, v = v, at = cbind(0.5, 0.5),
      vgm = list(model = "exponential", nugget = 0.01, psill = 0.3,
                 range = 0.4))
#>       pred        var
#> 1 1.062215 0.04918238
```
