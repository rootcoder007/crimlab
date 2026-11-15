# Sample points, by local density or uniformly

Draws `size` of the points without replacement. With `dens = TRUE` each
draw picks a point with probability proportional to its local density
`sum_j exp(-(d_ij / b)^2)`, `b` the median distance between points
(great-circle distance for EPSG:4326), so that dense areas are sampled
more; otherwise every point is equally likely.

## Usage

``` r
spsamp(x, size, dens = TRUE, seed = NULL)
```

## Arguments

- x:

  A
  [`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md)
  object.

- size:

  The number of points: one whole number from 0 to the number of points.

- dens:

  `TRUE` to sample by local density.

- seed:

  `NULL` (a seed drawn from R's generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes the draws) or
  one whole number from 0.

## Value

An integer vector of row numbers, with the attribute `prob`, the
probability of each point at the first draw.

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
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md),
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md),
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)

## Examples

``` r
s <- spobj(expand.grid(lon = 1:5, lat = 1:5))
spsamp(s, 5, seed = 1)
#> [1]  8 20 10  6 18
#> attr(,"prob")
#>  [1] 0.02813274 0.03613227 0.03913468 0.03613227 0.02813274 0.03614678
#>  [7] 0.04642113 0.05027751 0.04642113 0.03614678 0.03916871 0.05029579
#> [13] 0.05447252 0.05029579 0.03916871 0.03618351 0.04645522 0.05031121
#> [19] 0.04645522 0.03618351 0.02818869 0.03618420 0.03918600 0.03618420
#> [25] 0.02818869
```
