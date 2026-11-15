# Summarise a change table

One row per series (per group of a grouped table): the first and last
period with a value, the values there, the change from the first to the
last in percent, the compound rate per step, and how many of the
period-to-period changes went up, down or nowhere.

## Usage

``` r
yoysmry(x, lag = 1)
```

## Arguments

- x:

  A change table
  ([`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)).

- lag:

  The periods per step, the `lag` the table was made with.

## Value

A data frame of class `yoysmry`: `group` (for a grouped table), `from`,
`to`, `first`, `last`, `totalpct`, `cagrpct`, `periods`, `up`, `down`
and `flat`.

## Details

The series are put in period order first. `totalpct` and `cagrpct` are
`NA` unless the first value is positive (a compound rate out of zero
does not exist); `cagrpct` is \\100((last/first)^{1/s} - 1)\\ with \\s\\
= (values - 1) / `lag` steps.

## See also

[`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md),
[`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md).

## Examples

``` r
# simulated police-reported counts of two offences
d <- data.frame(year = rep(2019:2023, 2),
                offence = rep(c("Assault", "Robbery"), each = 5),
                n = c(812, 845, 790, 901, 954, 120, 131, 98, 105, 140))
yoysmry(yoypx(d$n, d$year, d$offence))
#> yoysmry: 2 series
#>    group from   to first last totalpct  cagrpct periods up down flat
#>  Assault 2019 2023   812  954 17.48768 4.111352       5  3    1    0
#>  Robbery 2019 2023   120  140 16.66667 3.928988       5  3    1    0
```
