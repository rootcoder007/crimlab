# Severity-weighted crime index

Police-reported crime weighted by seriousness, as the Crime Severity
Index of Statistics Canada is built (Wallace et al. 2009): each
offence's count is multiplied by its weight (there, the incarceration
rate of the offence times its mean prison sentence in days), the
products are summed for each unit (a city, a neighbourhood) and period,
divided by the population and expressed against a base cell as 100.

## Usage

``` r
crsev(
  x,
  off,
  w,
  unit = NULL,
  time = NULL,
  pop = NULL,
  per = 1e+05,
  base = NULL,
  sent = NULL
)
```

## Arguments

- x:

  The counts: numbers from 0, one per row.

- off:

  The offence of each row: labels.

- w:

  The weight of each offence: numbers from 0, named by offence (names of
  offences that do not occur are ignored) or in the order of the
  offences (a factor's levels, or the sorted labels).

- unit:

  `NULL` or the unit of each row: labels.

- time:

  `NULL` or the period of each row: labels.

- pop:

  `NULL` or the population of each row's unit and period, the same on
  all its rows: positive numbers.

- per:

  The rates are per `per` people: one positive number.

- base:

  `NULL`, or the row of the result whose index is 100.

- sent:

  `NULL`, or the mean sentence of each offence in the same form as `w`;
  the weights are then `w * sent` (an incarceration rate times a mean
  sentence length).

## Value

An object of class `crsev`, a list: `tab`, a data frame with one row per
unit and period (in the order of the units, then the periods): `unit`
and `time` when given, `count` (the counts summed), `sev` (the weighted
sum), with `pop` also `pop`, `rate` and `srate` (the count and the
weighted sum per `per` people), and with `base` also `idx` (`srate`, or
`sev` without `pop`, as 100 times its ratio to the base row's); `w`, the
weights used (`off`, `w`); `per` and `base`.

## Details

Give one row per unit, period and offence (rows that repeat a unit,
period and offence add up); an offence with no row in a cell counts 0
there. Beside the weighted rate the plain count rate shows whether the
mix of offences moved, not only their number.

## References

Wallace, M., Turner, J., Matarazzo, A. and Babyak, C. (2009). *Measuring
Crime in Canada: Introducing the Crime Severity Index and Improvements
to the Uniform Crime Reporting Survey*. Statistics Canada, catalogue
85-004-X.

## Examples

``` r
d <- data.frame(city = rep(c("A", "B"), each = 2),
                off = rep(c("theft", "robbery"), 2),
                n = c(10, 2, 8, 3), pop = 1000)
crsev(d$n, d$off, c(theft = 1, robbery = 50), unit = d$city,
      pop = d$pop, base = 1)
#> crsev: severity-weighted counts and rates per 1e+05, index 100 at row 1
#>  unit count sev  pop rate srate      idx
#>     A    12 110 1000 1200 11000 100.0000
#>     B    11 158 1000 1100 15800 143.6364
```
