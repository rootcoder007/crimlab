# Veil-of-darkness test

Compares each group's share of stops made in daylight, when a driver's
race can be seen, with its share in darkness, when it cannot. The ratio
of the two odds against a reference group is the factor by which
visibility changes that group's chance of being stopped; 1 means no
change.

## Usage

``` r
vodtst(day, night = NULL, ref = 1, shift = NULL)
```

## Arguments

- day:

  Stops by group in daylight: non-negative numbers, as a numeric vector,
  a 1-d table, a list of single numbers or a one-row data frame. With
  `night` left out, one table of groups by stops: a matrix, 2-d table or
  data frame whose two columns are named `day` and `night`, in either
  order (columns are never taken by position).

- night:

  Stops by group in darkness: positive numbers, in the forms of `day`;
  named inputs are matched to `day` by name, and an error lists any
  group named in one and not the other. Unnamed inputs are matched by
  position: only names protect against a mismatch.

- ref:

  The reference group, by name (matched exactly, case included; a factor
  by its label) or by position in `day` as given.

- shift:

  `NULL`, or one positive number per group: daylight exposure over
  darkness exposure, in the forms of `day`, matched by name.

## Value

A data frame with one row per group, in the order of `day`: `group`,
`ref` (`TRUE` for the reference group), `shday`, `shnight`, `ratio`,
`sign` (the sign of `shday - shnight`) and, when `shift` is given,
`adj`.

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
#>   group   ref shday shnight    ratio sign
#> 1 white  TRUE   0.6     0.7 1.000000   -1
#> 2 black FALSE   0.4     0.3 1.555556    1
# one table of groups by day and night
stops <- data.frame(day = c(600, 400), night = c(700, 300),
                    row.names = c("white", "black"))
vodtst(stops)
#>   group   ref shday shnight    ratio sign
#> 1 white  TRUE   0.6     0.7 1.000000   -1
#> 2 black FALSE   0.4     0.3 1.555556    1
```
