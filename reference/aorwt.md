# Aoristic weights

Spreads each event known only to have happened in an interval
`[from, to)` with uniform density over that interval and counts each
instant in the bin of its local clock time (Ratcliffe 2002): with the
defaults, the 168 hours of the week from Sunday 00:00. A change of clock
(daylight saving) moves the remaining time to the right hours: skipped
hours get nothing and a repeated hour gets both. Missing values follow
stated rules: an event with no `from` is left out; with no `to`, a
reversed interval or one shorter than `lo` seconds, its whole weight
goes to the start's bin; one of at least `hi` seconds is spread evenly
over all bins.

## Usage

``` r
aorwt(
  from,
  to,
  bw = 3600,
  nb = 168,
  origin = 259200,
  off = 0,
  tt = NULL,
  toff = NULL,
  lo = 120,
  hi = 604800,
  each = FALSE
)
```

## Arguments

- from, to:

  The start and end of each event's interval: date-times (POSIXct) or
  seconds since 1970-01-01 UTC; `NA` allowed as above, not infinite
  values.

- bw:

  The bin width in seconds, positive.

- nb:

  The number of bins in a period: one whole number from 1.

- origin:

  The local clock time at which bin 1 starts, in the same units (the
  default is Sunday 1970-01-04 00:00).

- off:

  The offset of the local clock from UTC in seconds, before `tt[1]` (for
  a fixed offset, the only one).

- tt, toff:

  `NULL`, or the instants at which the local clock's offset changes
  (strictly increasing) and the new offsets.

- lo, hi:

  The lengths in seconds below which an interval counts at its start and
  from which it counts everywhere, `0 <= lo < hi` (`lo = 0`, `hi = Inf`
  spread every interval exactly).

- each:

  `TRUE` to return each event's weights.

## Value

A list: `tot` (the summed weights by bin), `w` (`NULL`, or a matrix with
one row per event and one column per bin, rows summing to 1, 0 for
events left out) and `cnt` (events left out, without end, reversed,
short and long).

## References

Ratcliffe, J. H. (2002). Aoristic signatures and the spatio-temporal
analysis of high volume crime patterns. *Journal of Quantitative
Criminology*, 18(1), 23-43.
[doi:10.1023/A:1013240828824](https://doi.org/10.1023/A%3A1013240828824)

## Examples

``` r
f <- as.POSIXct(c("2024-03-04 22:00", "2024-03-05 09:30"), tz = "UTC")
t <- as.POSIXct(c("2024-03-05 02:00", "2024-03-05 10:30"), tz = "UTC")
a <- aorwt(f, t, bw = 3600, nb = 24, origin = 0)
a$tot
#>  [1] 0.25 0.25 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.50 0.50 0.00 0.00 0.00 0.00
#> [16] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.25 0.25
```
