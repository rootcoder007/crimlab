# CUSUM surveillance of counts

Watches a series of counts for a rise above what is expected: the
cumulative sum \\s_t = \max(0, s\_{t-1} + z_t - k)\\ of the standardised
excess \\z_t\\ raises an alarm where \\s_t \ge h\\. The expected counts
are given (`expd`), or are the mean of the counts before `start`, the
first period watched. `z` is `(y - m) / sqrt(m)` (`"std"`), Rossi's
transform `(y - 3 m + 2 sqrt(y m)) / (2 sqrt(m))` (`"rossi"`),
Anscombe's `1.5 (y^(2/3) - m^(2/3)) / m^(1/6)` (`"anscombe"`) or the
count itself (`"none"`, for `k` and `h` on the count scale). With
`reset` the sum starts again from 0 after an alarm.

## Usage

``` r
cusrv(
  y,
  expd = NULL,
  start = NULL,
  k = 1.04,
  h = 2.26,
  trans = "std",
  reset = FALSE
)
```

## Arguments

- y:

  The counts: numbers from 0.

- expd:

  `NULL`, or the expected count of each period: positive numbers.

- start:

  Without `expd`, the first period watched: one whole number from 2.

- k:

  The reference value subtracted each period: one number.

- h:

  The decision limit: one positive number.

- trans:

  The standardisation: `"std"`, `"rossi"`, `"anscombe"` or `"none"`.

- reset:

  `TRUE` to start again after each alarm.

## Value

A data frame, one row per period watched: `i`, its position; `y`;
`expd`; `z`; the sum `s`; and `alarm`.

## References

Page, E. S. (1954). Continuous inspection schemes. *Biometrika*,
41(1-2), 100-115.
[doi:10.1093/biomet/41.1-2.100](https://doi.org/10.1093/biomet/41.1-2.100)

Rossi, G., Lampugnani, L. and Marchi, M. (1999). An approximate CUSUM
procedure for surveillance of health events. *Statistics in Medicine*,
18(16), 2111-2122.
[doi:10.1002/(SICI)1097-0258(19990830)18:16\<2111::AID-SIM171\>3.0.CO;2-Q](https://doi.org/10.1002/%28SICI%291097-0258%2819990830%2918%3A16%3C2111%3A%3AAID-SIM171%3E3.0.CO%3B2-Q)

## Examples

``` r
y <- c(10, 12, 9, 11, 10, 13, 19, 22, 25, 12)
cusrv(y, start = 6)
#>    i  y expd         z        s alarm
#> 1  6 13 10.4 0.8062258 0.000000 FALSE
#> 2  7 19 10.4 2.6667468 1.626747 FALSE
#> 3  8 22 10.4 3.5970073 4.183754  TRUE
#> 4  9 25 10.4 4.5272678 7.671022  TRUE
#> 5 10 12 10.4 0.4961389 7.127161  TRUE
cusrv(y, expd = rep(11, 10), trans = "rossi", reset = TRUE)
#>     i  y expd          z        s alarm
#> 1   1 10   11 -0.3051028 0.000000 FALSE
#> 2   2 12   11  0.2982325 0.000000 FALSE
#> 3   3  9   11 -0.6181361 0.000000 FALSE
#> 4   4 11   11  0.0000000 0.000000 FALSE
#> 5   5 10   11 -0.3051028 0.000000 FALSE
#> 6   6 13   11  0.5904378 0.000000 FALSE
#> 7   7 19   11  2.2483195 1.208320 FALSE
#> 8   8 22   11  3.0321034 3.200423  TRUE
#> 9   9 25   11  3.7939546 2.753955  TRUE
#> 10 10 12   11  0.2982325 0.000000 FALSE
```
