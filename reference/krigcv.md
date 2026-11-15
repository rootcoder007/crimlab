# Leave-one-out cross-validation of ordinary kriging

Predicts each value from all the others by ordinary kriging with the
exponential covariance of
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
and summarises the residuals by the mean squared prediction error
`MSPE = (1 / n) sum (Z(s_i) - Z_-i(s_i))^2`, its root and the mean
absolute error (Schabenberger and Gotway 2005, ch. 4). Comparing it
across covariance parameters chooses among them.

## Usage

``` r
krigcv(x, y = NULL, v, nugget = 0, sill = 1, range = 1)
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

- nugget, sill, range:

  The covariance: `nugget` from 0, `sill` at least `nugget`, `range`
  positive.

## Value

A list of class `krigcv`: `mspe`, `rmspe`, `mae` and `resid` (the
leave-one-out residuals).

## References

Schabenberger, O. and Gotway, C. A. (2005). *Statistical Methods for
Spatial Data Analysis*. Chapman and Hall/CRC, Boca Raton.
[doi:10.1201/9781315275086](https://doi.org/10.1201/9781315275086)

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
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
set.seed(1)
p <- cbind(x = runif(30), y = runif(30))
krigcv(p, v = rnorm(30), range = 0.2)
#> krigcv: leave-one-out kriging of 30 values
#>      mspe     rmspe       mae 
#> 0.7067174 0.8406648 0.7098562 
```
