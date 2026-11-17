# Placement volatility events by fiscal year

The year-over-year table of placement volatility (Ruhela 2026, Appendix
E, Table E15) from the person-year table: for each fiscal year the
regional changes (`vm` summed), placements, person-years, the rates of
changes per 1,000 placements and per 1,000 person-years with exact
Poisson intervals, each year's share of all changes, and the percent
changes from the previous year of the count and of both rates. The rate
intervals are `qgamma(a/2, n)` to `qgamma(1 - a/2, n + 1)` over the
exposure, so a count of 0 has lower limit 0.

## Usage

``` r
cxyoy(x, level = 0.95, ...)
```

## Arguments

- x:

  A person-year table from
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md),
  or what
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  takes (a data frame of the release or a table key).

- level:

  The confidence level of the gamma intervals of both rates: one number
  from 0 to 1.

- ...:

  Passed to
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  when `x` is not a person-year table.

## Value

A data frame of class `cxyoy`, one row per fiscal year: `year`, `vm`,
`np`, `py`, `rnp`, `rnplo`, `rnphi` (per 1,000 placements), `rpy`,
`rpylo`, `rpyhi` (per 1,000 person-years), `share`, `dvm`, `drnp`,
`drpy` (percent changes from the previous year).

## Details

The two rates can disagree in direction: placements per person-year move
too, and the placement denominator absorbs their growth.

## References

Ruhela, V. S. (2026). Alert complexity and placement volatility in
Ontario restrictive confinement data. Major research paper, Centre for
Criminology and Sociolegal Studies, University of Toronto.

## Examples

``` r
p <- data.frame(year = rep(2024:2026, c(50, 60, 70)), vm = rpois(180, 0.2),
                np = rpois(180, 20) + 1)
cxyoy(p)
#> cxyoy: regional changes per fiscal year
#>  year vm   np py   rnp rnplo rnphi   rpy rpylo rpyhi  share     dvm   drnp
#>  2024 16 1064 50 15.04 8.595 24.42 320.0 182.9 519.7 0.3556      NA     NA
#>  2025 14 1218 60 11.49 6.284 19.29 233.3 127.6 391.5 0.3111 -12.500 -23.56
#>  2026 15 1493 70 10.05 5.623 16.57 214.3 119.9 353.4 0.3333   7.143 -12.59
#>     drpy
#>       NA
#>  -27.083
#>   -8.163
```
