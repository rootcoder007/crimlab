# The tailwind component of the wind along a flight course

The part of the wind that blows along a course, positive for a tailwind,
from which Hener (2022) builds the instrument of the landing direction
at Frankfurt Airport: aircraft land westbound (course 250 degrees)
unless the tailwind on that course exceeds 2.6 m/s, so the days of the
eastbound approach, and which municipalities lie under the noise, follow
from the wind.

## Usage

``` r
tailwd(dir, spd, course = 250, by = NULL, agg = "max")
```

## Arguments

- dir:

  Wind directions in degrees (`NA` allowed).

- spd:

  Wind speeds, one per direction, 0 or more (`NA` allowed).

- course:

  The course flown, degrees clockwise from north.

- by:

  `NULL` or the group of each reading.

- agg:

  `"max"` or `"mean"`: a group's value.

## Value

Without `by`, the tailwind component of each reading; with `by`, a data
frame, one row per group (sorted): `by` and `tw` (`NA` for a group
without a reading).

## Details

For wind from `dir` (degrees clockwise from north, where the wind blows
from) at speed `spd`, the component along `course` is
`spd * cos(dir + 180 - course)`: "positive for tailwinds and negative
for headwinds in the direction of the West25 landing approach", 0 for
"no wind at all or if it is blowing exactly perpendicular to the runway
(i.e., from 340 or 160 degrees)" (Hener 2022, sec. 3, for
`course = 250`). With `by` (the day of hourly readings), each group
gives one value: `agg = "max"`, the paper's "daily maximum partial wind
speed ... calculated from hourly wind data", or `"mean"`. The instrument
is then 1 for a municipality under the eastbound approach on a day whose
tailwind exceeds 2.6 m/s, 1 for one under the westbound approach on a
day whose tailwind is at most 2.6, and 0 otherwise (sec. 3; see
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md) for
the estimate).

## References

Hener, T. (2022). Noise pollution and violent crime. *Journal of Public
Economics*, 215, 104748.
[doi:10.1016/j.jpubeco.2022.104748](https://doi.org/10.1016/j.jpubeco.2022.104748)

## See also

[`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md),
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`dwnexp()`](https://rootcoder007.github.io/crimlab/reference/dwnexp.md)

## Examples

``` r
# wind from the east (70 degrees) at 4 m/s is a 4 m/s tailwind westbound
tailwd(c(70, 340, 250), c(4, 4, 4))
#> [1]  4.000000e+00 -7.347881e-16 -4.000000e+00
# hourly readings to the daily maximum, then Hener's instrument for a
# municipality under the eastbound approach
tw <- tailwd(c(70, 90, 250, 260), c(3, 2, 5, 1), by = c(1, 1, 2, 2))
as.integer(tw$tw > 2.6)
#> [1] 1 0
```
