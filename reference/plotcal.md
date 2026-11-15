# Many series over time: calendar, horizon, small-multiple and bump charts

`plotcal()` draws daily values as a calendar: one row of weeks per year,
each day a cell shaded by its value. `plothrz()` draws a horizon chart:
each series' distance from its own median, folded into `bands` bands of
darkening shade, blue above and red below, one row per series.
`plotmult()` draws each series in its own panel on a shared scale (small
multiples). `plotbump()` draws the rank of each series in each period (a
bump chart).

## Usage

``` r
plotcal(x, y, ..., main = NULL, col = NULL, pal = NULL)

plothrz(
  y,
  prd,
  grp = NULL,
  bands = 3L,
  ...,
  main = NULL,
  col = NULL,
  pal = NULL
)

plotmult(y, prd, grp = NULL, ..., main = NULL, col = NULL, pal = NULL)

plotbump(y, prd, grp, desc = TRUE, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x:

  For `plotcal()`, dates (or text
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html) reads), one per
  value of `y`.

- y:

  Numbers.

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

- prd:

  The period of each value: numbers, dates, a factor or text (in order
  of first appearance).

- grp:

  `NULL` (one series) or the series of each value; required by
  `plotbump()`.

- bands:

  The number of bands: one whole number from 1 to 6.

- desc:

  `TRUE` to rank the largest value first.

## Value

Invisibly, a data frame of what was drawn: for `plotcal()` the `date`,
`year`, `week`, `wday` (1 for Monday) and `y`; for the others `grp`,
`prd` and `y`, with `band` (`plothrz()`) or `rank` (`plotbump()`).

## Examples

``` r
days <- seq(as.Date("2024-01-01"), as.Date("2024-03-31"), by = "day")
plotcal(days, seq_along(days) %% 7)

d <- expand.grid(year = 2015:2020, area = c("n", "s", "e"))
d$v <- c(1:6, 6:1, rep(3, 6))
plothrz(d$v, d$year, d$area, bands = 2)

plotmult(d$v, d$year, d$area)

plotbump(d$v, d$year, d$area)
```
