# Small multiples of a choropleth map by period

`smmap()` draws one choropleth map per period in a grid of panels, all
with the same class breaks and colours, so that a change between periods
is a change of shade. The areas, breaks and classes are those of
[`chomap()`](https://rootcoder007.github.io/crimlab/reference/chomap.md),
the breaks computed from the values of every period together; the legend
is drawn in the last panel.

## Usage

``` r
smmap(
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

  A numeric matrix (or a data frame of numeric columns), one row per
  area and one column per period; the column names label the periods.

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

Invisibly, a list: `areas`, a data frame with one row per area and
period (the periods in order, the areas within each) holding `area`,
`prd` (the period's label), `v`, `cls` and `col`; `breaks`; and `bbox`,
as for
[`chomap()`](https://rootcoder007.github.io/crimlab/reference/chomap.md).

## Examples

``` r
sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))
p <- list(sq(0, 0), sq(1, 0), sq(0, 1), sq(1, 1))
v <- cbind(`2022` = c(3, 8, 5, 1), `2023` = c(4, 9, 2, 1),
           `2024` = c(6, 12, 2, 0))
smmap(p, v, brk = 3)
```
