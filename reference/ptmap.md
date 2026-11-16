# Point maps, with clusters

`ptmap()` draws points on a map at an aspect of 1, optionally coloured
by cluster: a result of
[`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md)
or a cluster number per point. Each cluster has its own colour and point
shape (so the clusters read in grey print too) and its number at its
centre, the mean of its points; noise is drawn as grey crosses. The box
and the cluster centres are computed by the compiled kernel; R only
draws. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
[`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md)
result draws it through `ptmap()`, given the points it was found in.

## Usage

``` r
ptmap(x, y = NULL, cl = NULL, ..., main = NULL, col = NULL, pal = NULL)

# S3 method for class 'dbscn'
plot(x, y = NULL, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x:

  The points: a numeric matrix or data frame of two columns (x and y),
  or the x coordinates with `y`. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), a result of
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md).

- y:

  `NULL` or the y coordinates, one per value of `x`. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the points
  the clusters were found in, as `x` of `ptmap()`.

- cl:

  `NULL`, a result of
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md),
  or the cluster of each point: a whole number from 1, `NA` for noise.

- ...:

  Passed to
  [`graphics::plot.window()`](https://rdrr.io/r/graphics/plot.window.html).

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours, one per cluster.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

Invisibly, a list: `pts`, a data frame of the points (`x`, `y`, `cl` the
cluster or `NA`, `col` and `pch` as drawn); `cent`, one row per cluster
number from 1 to the largest (`cl`, the centre `x` and `y`, `NaN` for a
number no point has, and `n` its points); and `bbox`, the box drawn
(`xmin`, `xmax`, `ymin`, `ymax`).

## Examples

``` r
set.seed(1)
p <- rbind(cbind(rnorm(30, 0, 0.2), rnorm(30, 0, 0.2)),
           cbind(rnorm(30, 3, 0.2), rnorm(30, 1, 0.2)),
           cbind(runif(10, -2, 5), runif(10, -2, 3)))
ptmap(p)

d <- dbscn(p, eps = 0.3, minpts = 5)
plot(d, p)
```
