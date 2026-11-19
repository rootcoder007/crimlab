# Spatial linear model by profile maximum likelihood

The Gaussian spatial linear model `v = Z beta + delta + e`, `delta` a
Gaussian process with covariance `sigma2 exp(-h / phi)`, without a
nugget (Schabenberger and Gotway 2005, ch. 5). The range `phi` is found
by maximising the profile likelihood with Brent's method on `log(phi)`
over `[log(max(h / 100, 0.001)),` `log(max(3 h, 1))]`, `h` the largest
distance, as R's [`optimize()`](https://rdrr.io/r/stats/optimize.html)
runs it; `beta` follows by generalised least squares. The standard
errors use `sigma2` estimated with `n - p` degrees of freedom.

## Usage

``` r
sglmx(v, z = NULL, x, y = NULL)
```

## Arguments

- v:

  The outcome: finite numbers, one per point.

- z:

  The design: `NULL` for an intercept alone, or a numeric matrix (or
  vector, or data frame) with one row per point; an intercept is not
  added.

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

## Value

A list of class `sglmx`: `coef`, `se`, `sigma2`, `phi`, `tau2` (the
nugget, 0) and `n`. [`summary()`](https://rdrr.io/r/base/summary.html)
gives the table of coefficients.

## References

Schabenberger, O. and Gotway, C. A. (2005). *Statistical Methods for
Spatial Data Analysis*. Chapman and Hall/CRC, Boca Raton.
[doi:10.1201/9781315275086](https://doi.org/10.1201/9781315275086)

## See also

Other geostatistics:
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md),
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md),
[`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md),
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
p <- cbind(x = runif(40), y = runif(40))
z <- cbind(1, p[, 1])
sglmx(1 + 2 * p[, 1] + rnorm(40, sd = 0.2), z, p)
#> sglmx: spatial linear model, exponential covariance, 40 points
#>         z1         z2     sigma2        phi 
#> 0.97349153 2.09618531 0.03478370 0.01495528 
```
