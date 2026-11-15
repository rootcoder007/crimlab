# Spatial linear mixed model by REML

The model `v = Z beta + delta + e`, `delta ~ N(0, sigma2 R_phi)` with
the exponential correlation `R_phi = exp(-h / phi)` and a nugget
`e ~ N(0, tau2 I)`, `tau2 = nu sigma2`, by restricted maximum likelihood
(Patterson and Thompson 1971; Schabenberger and Gotway 2005, ch. 5). The
search is R's Nelder-Mead
([`optim()`](https://rdrr.io/r/stats/optim.html), 400 evaluations,
relative tolerance `1e-6`) on `(log phi, log nu)` from
`(log(h / 3), log(0.1))`, `h` the largest distance (at least 1).

## Usage

``` r
smixm(v, z = NULL, x, y = NULL)
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

A list of class `smixm`: `coef`, `se`, `sigma2`, `tau2`, `phi`, `conv`
and `n`. [`summary()`](https://rdrr.io/r/base/summary.html) gives the
table of coefficients.

## References

Patterson, H. D. and Thompson, R. (1971). Recovery of inter-block
information when block sizes are unequal. *Biometrika*, 58(3), 545-554.
[doi:10.1093/biomet/58.3.545](https://doi.org/10.1093/biomet/58.3.545)

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
smixm(1 + 2 * p[, 1] + rnorm(40, sd = 0.3), cbind(1, p[, 1]), p)
#> smixm: spatial linear mixed model by REML, 40 points
#>         z1         z2     sigma2       tau2        phi 
#> 0.96148075 2.14380681 0.06465755 0.01544588 0.03247176 
```
