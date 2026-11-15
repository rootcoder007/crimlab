# Polynomial trend surface

Fits the mean of the values as a polynomial in the coordinates by
ordinary least squares, `mu(s) = sum_k beta_k f_k(s)`, the `f_k` the
monomials up to degree `order`: `1`; `x`, `y`; `x2`, `y2`, `xy`; `x3`,
`y3`, `x2y`, `xy2` (Schabenberger and Gotway 2005, ch. 2). The large-
scale trend is often removed this way before a variogram is fitted to
the residuals.

## Usage

``` r
trsurf(x, y = NULL, v, order = 2)
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

- order:

  The degree: 0, 1, 2 or 3.

## Value

A list of class `trsurf`: `coef` and `se` (named by monomial), `r2`,
`order` and `n`. [`summary()`](https://rdrr.io/r/base/summary.html)
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
[`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md),
[`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md),
[`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md),
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md),
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
v <- 1 + 2 * p[, 1] - p[, 2]^2 + rnorm(50, sd = 0.1)
f <- trsurf(p, v = v)
f
#> trsurf: trend surface of order 2 on 50 points, R-squared 0.9818
#>           1           x           y          x2          y2          xy 
#>  1.09625518  2.11305475 -0.73267272 -0.08022062 -0.33750250  0.14443429 
summary(f)
#>   term         est         se          z
#> 1    1  1.09625518 0.07968901 13.7566670
#> 2    x  2.11305475 0.21086371 10.0209505
#> 3    y -0.73267272 0.25326074 -2.8929581
#> 4   x2 -0.08022062 0.18258998 -0.4393484
#> 5   y2 -0.33750250 0.21097762 -1.5997076
#> 6   xy  0.14443429 0.18005051  0.8021877
```
