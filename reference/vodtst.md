# Veil-of-darkness test

Compares each group's share of stops made in daylight, when a driver's
race can be seen, with its share in darkness, when it cannot. The ratio
of the two odds against a reference group is the factor by which
visibility changes that group's chance of being stopped; 1 means no
change.

## Usage

``` r
vodtst(day, night, ref = 1, shift = NULL)
```

## Arguments

- day:

  Stops by group in daylight: non-negative numbers, optionally named.

- night:

  Stops by group in darkness: positive numbers, same length.

- ref:

  The reference group, by name (matched exactly, case included) or
  position.

- shift:

  `NULL`, or one positive number per group: daylight exposure over
  darkness exposure.

## Value

A data frame with one row per group: `group`, `share_day`,
`share_night`, `ratio`, `sign` (the sign of `share_day - share_night`)
and, when `shift` is given, `adj`.

## Details

The ratio identifies the visibility effect only when the mix of drivers
on the road is the same in daylight and darkness. `shift` states how
much each group's exposure differs between the two, and `adj` is the
ratio with that difference taken out. See `proof("vodtst")`.

## References

Grogger, J. and Ridgeway, G. (2006). Testing for racial profiling in
traffic stops from behind a veil of darkness. *Journal of the American
Statistical Association*, 101(475), 878-887.
[doi:10.1198/016214506000000168](https://doi.org/10.1198/016214506000000168)

## Examples

``` r
vodtst(day = c(white = 600, black = 400),
       night = c(white = 700, black = 300))
#>   group share_day share_night    ratio sign
#> 1 white       0.6         0.7 1.000000   -1
#> 2 black       0.4         0.3 1.555556    1
```
