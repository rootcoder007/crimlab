# Outcome by adoption cohort

The mean outcome of each adoption cohort by time since its adoption,
smoothed by a moving average and scaled to 1 in the period before
adoption: a look, before any model, at whether the cohorts moved
together before the policy. The period before adoption is shaded.

## Usage

``` r
plotcoh(
  y,
  unit,
  time,
  first,
  win = c(-12, 12),
  ma = 3,
  norm = TRUE,
  ...,
  main = NULL,
  col = NULL,
  pal = NULL
)
```

## Arguments

- y:

  The outcome: numbers, one per row.

- unit, time:

  The unit and the period of each row; `time` numbers.

- first:

  The period the unit of each row is first treated, `NA` for a unit
  never treated (those are not drawn).

- win:

  The periods before and after adoption to show: two whole numbers, the
  first -1 or below and the second 0 or above.

- ma:

  The periods in the moving average: a whole number from 1 (no
  smoothing) to 100.

- norm:

  `TRUE` to divide each cohort by its value in the period before
  adoption.

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

## Value

Invisibly, a data frame of the points drawn: `cohort`, `evtm` (periods
since adoption) and `value`.

## Details

A cohort is the units first treated in the same period. Its mean in a
period is over its units observed then; the moving average is over that
period and the `ma - 1` before it, so the first `ma - 1` periods of a
cohort have no value. With `norm`, a cohort whose period before adoption
is not in the window, or is 0, has no values.

## Examples

``` r
d <- expand.grid(year = 2012:2024, town = LETTERS[1:12])
i <- match(d$town, LETTERS)
d$first <- rep(c(2018, 2021, NA), each = 4)[i]
d$rate <- 30 + i + 0.5 * (d$year - 2012) + sin(seq_len(156)) -
  4 * (!is.na(d$first) & d$year >= d$first)
plotcoh(d$rate, d$town, d$year, d$first, win = c(-4, 3), ma = 1)
```
