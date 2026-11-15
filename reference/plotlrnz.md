# Lorenz and funnel plots

`plotlrnz()` draws the Lorenz curve of a total spread over units: the
cumulative share of the total held by the cumulative share of units,
smallest first, against the line of equality, with the Gini coefficient
(from
[`cconc()`](https://rootcoder007.github.io/crimlab/reference/cconc.md))
in the title. `plotfnl()` draws a funnel plot (Spiegelhalter 2005): each
area's ratio of observed to expected counts against its expected count,
with exact Poisson control limits around `target`; areas outside the
outer limits are drawn in a second colour and a second shape, and
labelled when `area` is given.

## Usage

``` r
plotlrnz(x, ..., main = NULL, col = NULL, pal = NULL)

plotfnl(
  obs,
  expd,
  area = NULL,
  lvls = c(0.95, 0.998),
  target = 1,
  ...,
  main = NULL,
  col = NULL,
  pal = NULL
)
```

## Arguments

- x:

  For `plotlrnz()`, non-negative values, one per unit, or a data frame
  with columns `pop` and `share` (a Lorenz curve already computed, such
  as `lorenz` of
  [`cconc()`](https://rootcoder007.github.io/crimlab/reference/cconc.md)).

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

- obs, expd:

  Observed and expected counts, one per area; `expd` positive.

- area:

  `NULL` or a label per area.

- lvls:

  The levels of the control limits: probabilities strictly between 0 and
  1.

- target:

  The ratio expected under no excess: one positive number.

## Value

Invisibly, the data drawn: for `plotlrnz()` the curve (`pop`, `share`)
from the origin, with the Gini as attribute `gini` (`NA` for a curve
given as a frame); for `plotfnl()` a data frame with `area`, `observed`,
`expected`, `ratio` and `outside` (beyond the outer limits).

## Details

The funnel limits at expected count \\E\\ and level \\p\\ are the
Poisson quantiles \\q\_{(1-p)/2}\\ and \\q\_{(1+p)/2}\\ of a count with
mean `target` \\E\\, divided by \\E\\; they are drawn through the areas'
own expected counts. An area with either count missing is dropped whole.

## References

Spiegelhalter, D. J. (2005). Funnel plots for comparing institutional
performance. *Statistics in Medicine*, 24(8), 1185-1202.
[doi:10.1002/sim.1970](https://doi.org/10.1002/sim.1970)

## Examples

``` r
plotlrnz(c(0, 0, 1, 1, 2, 3, 8, 25))

plotfnl(obs = c(4, 20, 50, 2), expd = c(5, 22, 40, 20),
        area = c("a", "b", "c", "d"))
```
