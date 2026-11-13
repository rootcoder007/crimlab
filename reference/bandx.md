# Banded categories: bounds, values and what the open band costs

Justice statistics are often published in bands ("1-5 days", "16 to 30",
"Over 60", "Less than 1 year"). `bandx()` reads each label into a lower
and an upper end, places a value in each band (the midpoint, an end, or
the geometric mean), expands a banded frequency table into one value per
unit, and shows how much a statistic of those values depends on the cap
assumed for an open top band.

## Usage

``` r
bandx(
  x,
  cnts = NULL,
  rule = "midpoint",
  cap = NULL,
  fac = 2,
  floor = NULL,
  clsup = TRUE,
  intsca = TRUE,
  caps = NULL,
  stat = NULL,
  drop = FALSE
)
```

## Arguments

- x:

  The band labels: text or a factor.

- cnts:

  `NULL`, or the number of units in each band: numbers from 0.

- rule:

  Where in a band its value sits: `"midpoint"`, `"lower"`, `"upper"` or
  `"geometric"` (the midpoint where the lower end is 0).

- cap:

  `NULL`, or the end assumed for an open top band.

- fac:

  Without `cap`, the open top band ends at `fac` times its start: one
  positive number.

- floor:

  `NULL`, or the start assumed for an open bottom band.

- clsup:

  `TRUE` when `"60+"` includes 60.

- intsca:

  `TRUE` when the bands hold whole numbers, so that `"over 60"` starts
  at 61 and `"under 10"` ends at 9.

- caps:

  `NULL`, or the caps at which to recompute `stat` (by default 12 caps
  from 1.05 to 20 times the open band's start).

- stat:

  `NULL`, or a function of the expanded values giving one number (with
  `cnts`, and an open top band).

- drop:

  `TRUE` to leave out the units of bands that cannot be read; otherwise
  they are an error.

## Value

An object of class `bandx`, a list: `bands`, a data frame (`label`,
`lower`, `upper`, `openlo`, `openup`, `value`, `lo`, `hi` the ends
used); `vals`, the expanded values (with `cnts`); `sens` (with `stat`):
a data frame of `cap` and `value`, and `span` and `relspan` (the span
over the caps, absolute and over the median absolute value).

## Details

A label is a range (`"5-10"`, `"5 to 10"`, `"5 through 10"`; any dash),
a single number, or a number with a word that opens it: `"60+"`,
`"60 and over"`, `"60 or more"` (from 60), `"over 60"`, `"more than 60"`
(from 61 for whole numbers, `intsca`), `"10 or less"` (to 10),
`"under 10"`, `"less than 10"`, `"<10"` (to 9 for whole numbers).
Thousands separators are dropped. With `clsup = FALSE` a `"60+"` band
starts at 61. An open top band ends at `cap`, or at `fac` times its
start; an open bottom band starts at `floor`, or at 0.

## Examples

``` r
b <- c("1-5", "6-15", "16-30", "31-60", "Over 60")
n <- c(40, 25, 15, 12, 8)
r <- bandx(b, n, stat = mean)
r$bands
#>     label lower upper openlo openup value lo  hi
#> 1     1-5     1     5  FALSE  FALSE   3.0  1   5
#> 2    6-15     6    15  FALSE  FALSE  10.5  6  15
#> 3   16-30    16    30  FALSE  FALSE  23.0 16  30
#> 4   31-60    31    60  FALSE  FALSE  45.5 31  60
#> 5 Over 60    61    NA  FALSE   TRUE  91.5 61 122
r$sens
#>           cap    value
#> 1    64.05000 17.73700
#> 2    83.72753 18.52410
#> 3   109.45042 19.55302
#> 4   143.07594 20.89804
#> 5   187.03193 22.65628
#> 6   244.49214 24.95469
#> 7   319.60536 27.95921
#> 8   417.79498 31.88680
#> 9   546.15054 37.02102
#> 10  713.93968 43.73259
#> 11  933.27724 52.50609
#> 12 1220.00000 63.97500
```
