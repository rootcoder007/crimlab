# Distributions by group: ridgelines, distribution functions, violins, hexagonal bins and concentration curves

Plots of many values, split by a group. `plotrdg()` stacks a kernel
density per group (a ridgeline plot); `plotcdfs()` overlays the
empirical distribution functions; `plotvln()` draws a violin per group
with its median and quartiles; `plothex()` counts two variables in
hexagonal bins, shaded by count; `plotgain()` overlays the concentration
(Lorenz) curves of the groups with their Gini coefficients.

## Usage

``` r
plotrdg(x, grp = NULL, bw = "nrd0", ..., main = NULL, col = NULL, pal = NULL)

plotcdfs(x, grp = NULL, ..., main = NULL, col = NULL, pal = NULL)

plotvln(x, grp = NULL, bw = "nrd0", ..., main = NULL, col = NULL, pal = NULL)

plothex(x, y, bins = 30L, ..., main = NULL, col = NULL, pal = NULL)

plotgain(x, grp = NULL, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x:

  Numbers (non-negative for `plotgain()`).

- grp:

  `NULL` (one group) or a group per value of `x`.

- bw:

  The kernel bandwidth: one positive number or a rule known to
  [`stats::density()`](https://rdrr.io/r/stats/density.html) (`"nrd0"`,
  `"nrd"`, `"ucv"`, `"bcv"`, `"sj"`, `"sj-ste"`, `"sj-dpi"`).

- ...:

  Passed to the base-graphics call that draws the frame.

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

- y:

  For `plothex()`, numbers, one per value of `x`.

- bins:

  The number of hexagons across: one whole number from 2 to 500.

## Value

Invisibly, a data frame of what was drawn: `grp`, `x` and `y` (the
density, or the share at or below), for `plotvln()` with the quartiles
as attribute `summary`; for `plothex()` the hexagon centres `x`, `y` and
`count`, by row of hexagons and then column; for `plotgain()` `grp`,
`pop` and `share`, with the Gini per group as attribute `gini`. Values
that are not finite are left out.

## Examples

``` r
set.seed(1)
x <- c(rnorm(100), rnorm(100, 2))
g <- rep(c("a", "b"), each = 100)
plotrdg(x, g)

plotcdfs(x, g)

plotvln(x, g)

plothex(x, x + rnorm(200), bins = 10)

plotgain(rexp(200), g)
```
