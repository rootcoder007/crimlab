# Constant reporters of a reporting panel

The units (agencies reporting to a national crime statistics programme,
say) that report at least `minrep` months in every year from `y0` to
`y1`. A year without a record of the unit, or with too few months,
disqualifies it, whatever other units do that year. Counts over time are
comparable only on such a constant panel; a trend from all reporters
mixes change in crime with change in who reports.

## Usage

``` r
rptc(unit, year, rep, y0, y1, minrep = 12)
```

## Arguments

- unit:

  The unit of each record: labels, no missing values.

- year:

  The year of each record: whole numbers.

- rep:

  The months reported in each record: whole numbers.

- y0, y1:

  The first and last years, `y0 <= y1`, from 0.

- minrep:

  The fewest months that count as reporting: a whole number from 0.

## Value

A data frame with one row per unit (sorted): `unit` and `const` (`TRUE`
for a constant reporter).

## References

Maltz, M. D. (1999). *Bridging Gaps in Police Crime Data*. Bureau of
Justice Statistics, NCJ 176365.

## Examples

``` r
rptc(unit = c("a", "a", "b", "b", "b"), year = c(2020, 2021, 2020, 2021,
     2021), rep = c(12, 12, 12, 6, 12), y0 = 2020, y1 = 2021)
#>   unit const
#> 1    a  TRUE
#> 2    b  TRUE
```
