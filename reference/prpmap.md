# Proportional symbol maps

`prpmap()` draws a circle at each point (a district's centre, say) whose
area is proportional to its value: the radius is `rmax` times the square
root of the value over the largest value, so the largest has radius
`rmax`. The circles are drawn largest first, so that the smaller stay
visible over them, and filled half-transparent. The radii, the drawing
order and the box (widened by `rmax` on every side so the circles fit)
are computed by the compiled kernel; R only draws.

## Usage

``` r
prpmap(x, y = NULL, v, rmax = NULL, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x:

  The points: a numeric matrix or data frame of two columns (x and y),
  or the x coordinates with `y`.

- y:

  `NULL` or the y coordinates, one per value of `x`.

- v:

  Numbers from 0, one per point; `NA` draws no circle (not all missing).

- rmax:

  `NULL` (a fifteenth of the longer side of the points' box) or the
  radius of the largest circle, in the units of the coordinates.

- ...:

  Passed to
  [`graphics::plot.window()`](https://rdrr.io/r/graphics/plot.window.html).

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or a colour.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

Invisibly, a list: `pts`, a data frame of the points (`x`, `y`, `v`, `r`
the radius, `NA` when not drawn, and `ord` the place in the drawing
order, `NA` when not drawn); `rmax`; and `bbox`, the box drawn (`xmin`,
`xmax`, `ymin`, `ymax`).

## Examples

``` r
x <- c(0, 4, 1, 3, 2)
y <- c(0, 1, 3, 4, 2)
prpmap(x, y, c(120, 30, 5, 60, 0))

prpmap(cbind(x, y), v = c(120, 30, NA, 60, 10), rmax = 1)
```
