# Case flow through a court

Cases received and disposed by period: the clearance rate (disposed over
received, above 1 when the caseload shrinks) with its exact interval,
the net, the clearance over a rolling window, and with the pending
caseload its change and the periods it would take to clear at the recent
pace. With the length of each period (`span`) the filing rate gives, by
Little's law, the mean time to disposition the pending caseload implies,
to set beside a published mean (`mtime`): their `ratio` far from 1 means
the published mean leaves out the long-pending cases. Published duration
bands (`shares` of cases in bands starting at `floors`) bound the mean
from below.

## Usage

``` r
dktflw(
  recv = NULL,
  disp = NULL,
  pend = NULL,
  prd = NULL,
  win = 12,
  level = 0.95,
  span = NULL,
  mtime = NULL,
  shares = NULL,
  floors = NULL,
  mix = NULL,
  scores = NULL
)
```

## Arguments

- recv, disp:

  `NULL`, or the cases received and disposed in each period: numbers
  from 0.

- pend:

  `NULL`, or the cases pending at the end of each period.

- prd:

  `NULL` (1, 2, ...) or the label of each period.

- win:

  The rolling window in periods: a whole number from 1.

- level:

  The level of the intervals: one number between 0 and 1.

- span:

  `NULL`, or the days in each period (one number or one per period);
  needs `pend`.

- mtime:

  `NULL`, or the published mean days to disposition of each period;
  needs `span`.

- shares, floors:

  `NULL`, or the share of cases in each duration band (a matrix, one row
  per period) and each band's lower edge.

- mix:

  `NULL`, or the outcome counts: a numeric matrix of periods by outcomes
  (named by its dimnames).

- scores:

  `NULL` (1, 2, ...) or the period scores of the trend test.

## Value

An object of class `dktflw`, a list: `flow` (with `recv`: `period`,
`recv`, `disp`, `clear`, `lo`, `hi`, `net`, `rclear`, `rlo`, `rhi`, then
`pend`, `pchg`, `toclear`, or `cumnet`, and with `span` `frate`, `impl`,
with `mtime` `ratio`, `ltlbk`, with `shares` `floorbd`); `mix` (a list
of `shares`, `homog` and `trend`); and `level`.

## Details

With `mix`, the counts of each outcome (conviction, withdrawal,
acquittal) by period: each period's shares with Goodman's simultaneous
intervals, the chi-square test that the mix is the same in every period
with Cramer's V, and the Cochran-Armitage trend in each share over the
periods (as
[`stats::prop.trend.test()`](https://rdrr.io/r/stats/prop.trend.test.html))
with Holm's adjustment.

## References

Goodman, L. A. (1965). On simultaneous confidence intervals for
multinomial proportions. *Technometrics*, 7(2), 247-254.
[doi:10.1080/00401706.1965.10490252](https://doi.org/10.1080/00401706.1965.10490252)

## Examples

``` r
dktflw(recv = c(120, 130, 125, 140), disp = c(110, 135, 120, 150),
       pend = c(400, 395, 400, 390), win = 2, span = 91)
#> dktflw: case flow
#>  period recv disp     clear        lo       hi net   rclear       rlo      rhi
#>       1  120  110 0.9166667 0.7011985 1.197368  10       NA        NA       NA
#>       2  130  135 1.0384615 0.8101053 1.331593  -5 0.980000 0.8183090 1.173541
#>       3  125  120 0.9600000 0.7409969 1.243285   5 1.000000 0.8373357 1.194264
#>       4  140  150 1.0714286 0.8453024 1.358733 -10 1.018868 0.8568263 1.211642
#>  pend pchg  toclear    frate     impl
#>   400   NA 3.636364 1.318681 303.3333
#>   395   -5 3.224490 1.428571 276.5000
#>   400    5 3.137255 1.373626 291.2000
#>   390  -10 2.888889 1.538462 253.5000
m <- matrix(c(50, 55, 60, 30, 28, 25, 20, 17, 15), 3,
            dimnames = list(2021:2023, c("guilty", "withdrawn", "other")))
dktflw(mix = m)$mix$trend
#>     outcome      stat         p      holm change
#> 1    guilty 2.0202020 0.1552185 0.4656555   0.10
#> 2 withdrawn 0.6246183 0.4293363 0.7006036  -0.05
#> 3     other 0.8723635 0.3503018 0.7006036  -0.05
```
