# Fit a variogram model by Gaussian maximum likelihood

Fits the nugget, partial sill and range of a spherical, exponential or
Gaussian covariance model by maximising the Gaussian log-likelihood of
the values about a constant mean (Mardia and Marshall 1984). The search
is the Nelder-Mead simplex of R's
[`optim()`](https://rdrr.io/r/stats/optim.html) (at most 800
evaluations) on the logarithms of the three parameters, from 0.1 and 0.9
of the variance of the values and a quarter of the largest distance. The
covariance matrix gets `1e-8` of the variance on its diagonal so that it
stays positive definite.

## Usage

``` r
vgfit(x, y = NULL, v, model = "exponential")
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

- model:

  `"exponential"`, `"spherical"` or `"gaussian"`.

## Value

A list of class `vgfit`: `model`, `nugget`, `psill`, `range`, `loglik`
(the log-likelihood less `n / 2 log(2 pi)`), `conv` (`TRUE` when the
simplex converged) and `nfn` (the evaluations).

## References

Mardia, K. V. and Marshall, R. J. (1984). Maximum likelihood estimation
of models for residual covariance in spatial regression. *Biometrika*,
71(1), 135-146.
[doi:10.1093/biomet/71.1.135](https://doi.org/10.1093/biomet/71.1.135)

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
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(40, 0, 10), y = runif(40, 0, 10))
v <- p[, 1] * 0.5 + rnorm(40, 0, 0.3)
vgfit(p, v = v)
#> vgfit: exponential model by maximum likelihood
#>       nugget        psill        range       loglik 
#>  0.006099213  8.719291632 71.862753310 11.408010219 
```
