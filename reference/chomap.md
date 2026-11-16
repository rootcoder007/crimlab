# Choropleth maps in base graphics

`chomap()` shades each area of a map by the class of its value, with no
spatial package: an area is a matrix of coordinates. The class breaks,
the class of each value and the bounding box are computed by the
compiled kernel; R only draws. Areas are drawn with
[`polypath()`](https://rdrr.io/r/graphics/polypath.html), so a ring
inside another is a hole.

## Usage

``` r
chomap(
  polys,
  v,
  brk = 5,
  meth = "quantile",
  ...,
  main = NULL,
  col = NULL,
  pal = NULL
)
```

## Arguments

- polys:

  A list of areas, each a numeric matrix of two columns (x and y) with
  at least three points; the rings of one area are separated by a row of
  `NA`. Names, when given, label the areas.

- v:

  Numbers, one per area.

- brk:

  The number of classes (a whole number from 1 to 9), or the breaks
  themselves (increasing numbers covering the finite values of `v`).

- meth:

  How a number of classes is cut: `"quantile"` or `"equal"`.

- ...:

  Passed to
  [`graphics::plot.window()`](https://rdrr.io/r/graphics/plot.window.html).

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours, recycled over the classes.

- pal:

  `NULL` (`"sequential"`) or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

Invisibly, a list: `areas`, a data frame with `area` (the label or the
position of each area), `v`, `cls` (its class from 1, `NA` when
unshaded) and `col` (its fill, `NA` when unshaded); `breaks`; and
`bbox`, the box drawn (`xmin`, `xmax`, `ymin`, `ymax`).

## Details

With `brk` a number of classes, `meth = "quantile"` puts the breaks at
R's type-7 quantiles of the finite values (so the classes hold about
equal counts; breaks that repeat are merged, which leaves fewer classes)
and `meth = "equal"` cuts their range into classes of equal width, as
`seq(min, max, length.out = brk + 1)`. A value falls in the class whose
interval holds it, intervals closed on the left and the last closed on
both sides, as
`findInterval(v, breaks, rightmost.closed = TRUE, all.inside = TRUE)`
gives. When every finite value is the same, the one class runs from 0.5
below it to 0.5 above. Values that are missing or not finite are left
unshaded.

The colours are a sequential ramp (light for low classes, dark for
high), so the order of the classes reads in grey print too.

## Examples

``` r
sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))
p <- list(a = sq(0, 0), b = sq(1, 0), c = sq(0, 1), d = sq(1, 1))
chomap(p, c(12, 3, 7, 25), brk = 2)

chomap(p, c(12, 3, 7, 25), brk = c(0, 5, 10, 30), main = "Rates")
```
