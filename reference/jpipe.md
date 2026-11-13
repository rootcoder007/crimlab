# Rates between stages of the justice pipeline

Cases moving from one stage to the next (charges to convictions,
admissions to releases) take time, so a rate of this period's outflow
over this period's inflow mixes cohorts. `jpipe()` divides the numerator
summed over a window by the denominator summed over the same window
moved back by the lag between the stages.

## Usage

``` r
jpipe(num, den, win = 12, lag = 0, grp = NULL)
```

## Arguments

- num:

  The counts of the later stage, in time order: numbers from 0.

- den:

  The counts of the earlier stage: numbers from 0, one per `num`.

- win:

  The window, in periods: a whole number from 1.

- lag:

  The lag between the stages, in periods: a whole number from 0.

- grp:

  `NULL`, or the group of each period (each group in time order):
  labels.

## Value

A data frame of class `jpipe` with one row per period, the groups in
label order: `group`, `index` (the period within its group), `num` and
`den` (the window sums, `NA` until a full window), `rate`.

## Examples

``` r
jpipe(num = c(5, 7, 6, 9, 8, 10), den = c(20, 22, 19, 25, 24, 26),
      win = 3, lag = 1)
#> jpipe: lagged rates between stages
#>  group index num den      rate
#>    all     1  NA  NA        NA
#>    all     2  NA  NA        NA
#>    all     3  18  NA        NA
#>    all     4  22  61 0.3606557
#>    all     5  23  66 0.3484848
#>    all     6  27  68 0.3970588
```
