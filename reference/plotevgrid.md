# Event studies for several outcomes in one grid

Draws the event-study path of each result in a named list, one panel
each, as [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
one: estimates by time since adoption with their intervals, before
adoption in open circles and after in filled ones.

## Usage

``` r
plotevgrid(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x:

  A named list of results of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  or
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md).

- ...:

  Passed to the base-graphics call that draws each panel;
  `style = "stata"` draws hollow points joined by a line with dashed
  intervals, in place of open points before adoption and filled ones
  after.

- main:

  `NULL` or one string, an overall title.

- col, pal:

  As for
  [`plotlrnz()`](https://rootcoder007.github.io/crimlab/reference/plotlrnz.md).

## Value

Invisibly, a data frame of the points drawn: `panel`, `t`, `est`, `lo`
and `hi`.

## Examples

``` r
d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
i <- match(d$town, LETTERS)
first <- rep(c(2019, 2021, NA), each = 4)
d$first <- first[i]
d$rate <- 30 + i + 0.5 * (d$year - 2015) -
  4 * (!is.na(d$first) & d$year >= d$first) + sin(seq_len(120))
d$cases <- d$rate / 2 + cos(seq_len(120))
plotevgrid(list(
  Rate = evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3),
  Cases = evtst(d$cases, d$town, d$year, d$first, leads = 3, lags = 3)
))
```
