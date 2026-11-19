# Ordinary block kriging

Predicts the mean value over each block (an area such as a district) by
ordinary kriging with point-to-block covariances averaged over the
block, `C(s_i, B) = (1 / |B|) int_B` `C(s_i, u) du`, under the
exponential covariance `(sill - nugget) exp(-h / range)`, plus `nugget`
at distance 0 (Schabenberger and Gotway 2005, ch. 4). A block is a box,
given as a 2 x 2 matrix of its lower and upper corners (rows) and
integrated on a regular grid of `round(sqrt(nq))^2` points, or a matrix
of the points that discretise it.

## Usage

``` r
bkrig(x, y = NULL, v, blocks, nq = 25, nugget = 0, sill = 1, range = 1)
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

- blocks:

  A list of blocks: each a 2 x 2 matrix (rows: lower and upper corner,
  columns: x and y), or a matrix of two columns of points.

- nq:

  The size of the grid of a box: one whole number from 1. A box is
  averaged over `round(sqrt(nq))^2` points (`nq` rounded to a square: 10
  gives 9). The points of all the blocks must be at most 67108864
  (2^26).

- nugget, sill, range:

  The covariance: `nugget` from 0, `sill` at least `nugget`, `range`
  positive.

## Value

A data frame with one row per block: `est` (the kriged mean) and `se`
(its kriging standard error).

## References

Schabenberger, O. and Gotway, C. A. (2005). *Statistical Methods for
Spatial Data Analysis*. Chapman and Hall/CRC, Boca Raton.
[doi:10.1201/9781315275086](https://doi.org/10.1201/9781315275086)

## See also

Other geostatistics:
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
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(30), y = runif(30))
v <- p[, 1] + rnorm(30, sd = 0.1)
bkrig(p, v = v, blocks = list(rbind(c(0, 0), c(0.5, 0.5)),
                              rbind(c(0.5, 0.5), c(1, 1))),
      range = 0.3)
#>         est        se
#> 1 0.3357770 0.1852905
#> 2 0.7159535 0.2114923
```
