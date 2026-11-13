# Year-over-year change in a published series

The change of each period against the period `lag` before it, within
groups (an agency, a region): the difference, the percent change (or,
for percentages, the change in points), and for counts the exact
interval of the percent change, from the conditional binomial of the
count now among both, the exact two-sided p-value of no change (as
[`stats::binom.test()`](https://rdrr.io/r/stats/binom.test.html)), and
the p-values adjusted for the number of comparisons (as
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)). With
`rnd`, the counts are taken as rounded, and the envelope of the change
they allow is added, with a percent interval that covers both.

## Usage

``` r
yoypx(
  x,
  prd,
  by = NULL,
  lag = 1,
  units = "count",
  minbase = 0,
  level = 0.95,
  rnd = NULL,
  kind = "nearest",
  adj = "BH",
  alpha = 0.05
)
```

## Arguments

- x:

  The values: counts from 0, or numbers (`units` other than `"count"`).

- prd:

  The period of each row: numbers (a year, a month index); the earlier
  period is `prd - lag`.

- by:

  `NULL`, or the group of each row: labels.

- lag:

  How many periods back: one positive number.

- units:

  `"count"`, `"continuous"` or `"percent"`.

- minbase:

  Changes from an earlier value below this are flagged: one number from
  0.

- level:

  The level of the intervals: one number between 0 and 1.

- rnd:

  `NULL`, or the rounding base of the counts.

- kind:

  `"nearest"` or `"random"` rounding.

- adj:

  The adjustment of the p-values: `"none"`, `"bonferroni"`, `"holm"`,
  `"hochberg"`, `"BH"` or `"BY"`.

- alpha:

  The level at which `signif` marks an adjusted p-value.

## Value

A data frame of class `yoypx` with one row per group and period (in that
order): `group` (with `by`), `prd`, `value`, `prev`, `change`, `pctchg`
(`ppchg` for percentages), `flag`; for counts also `pctlo`, `pctup`,
`p`, `padj`, `signif`, and with `rnd` `envlo`, `envup`, `envplo`,
`envpup` (the envelope of the change and of the percent change) and
`cpctlo`, `cpctup`.

## Details

Rows that share a group and a period are summed (counts) or averaged. A
flag marks a change with no earlier period, an earlier period of 0, or
one below `minbase`.

## Examples

``` r
d <- data.frame(year = rep(2019:2023, 2),
                svc = rep(c("A", "B"), each = 5),
                n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33))
yoypx(d$n, d$year, d$svc, rnd = 5)
#> yoypx: change on the earlier period
#>  group  prd value prev change     pctchg              flag      pctlo
#>      A 2019   120   NA     NA         NA no earlier period         NA
#>      A 2020   131  120     11   9.166667              <NA> -15.436293
#>      A 2021    98  131    -33 -25.190840              <NA> -43.017280
#>      A 2022   105   98      7   7.142857              <NA> -19.421106
#>      A 2023   140  105     35  33.333333              <NA>   2.787153
#>      B 2019    30   NA     NA         NA no earlier period         NA
#>      B 2020    28   30     -2  -6.666667              <NA> -46.269451
#>      B 2021    41   28     13  46.428571              <NA> -11.610446
#>      B 2022    25   41    -16 -39.024390              <NA> -64.464747
#>      B 2023    33   25      8  32.000000              <NA> -23.859990
#>       pctup          p      padj signif envlo envup     envplo    envpup
#>          NA         NA        NA     NA    NA    NA         NA        NA
#>   41.033275 0.52799664 0.7039955  FALSE     6    16   4.897959  13.61702
#>   -2.059398 0.03423381 0.1369352  FALSE   -38   -28 -28.464419 -21.78988
#>   42.569190 0.67377407 0.7700275  FALSE     2    12   1.990050  12.56545
#>   73.399950 0.02963430 0.1369352  FALSE    30    40  27.906977  39.02439
#>          NA         NA        NA     NA    NA    NA         NA        NA
#>   61.661111 0.89568321 0.8956832  FALSE    -7     3 -21.538462  10.90909
#>  145.863553 0.14803218 0.2960644  FALSE     8    18  26.229508  70.58824
#>    2.692906 0.06401750 0.1707133  FALSE   -21   -11 -48.275862 -28.57143
#>  131.538008 0.35814330 0.5730293  FALSE     3    13  10.909091  57.77778
#>      cpctlo     cpctup
#>          NA         NA
#>  -15.436293  41.033275
#>  -43.017280  -2.059398
#>  -19.421106  42.569190
#>    2.787153  73.399950
#>          NA         NA
#>  -46.269451  61.661111
#>  -11.610446 145.863553
#>  -64.464747   2.692906
#>  -23.859990 131.538008
```
