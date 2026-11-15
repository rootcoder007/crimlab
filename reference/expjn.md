# Join an exposure table to counts by place and day

Makes the analysis table of the estimators of this family: the rows of
the outcome table (place, day, count, ...) with the exposure of the same
place and day joined on. Exposure may come from monitors (assigned to
areas beforehand, for example by inverse squared distance from the three
nearest, as Bondy, Roth and Sager do) or from crimlab's dispersion
functions
([`gplume()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
[`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md)).

## Usage

``` r
expjn(expo, outc, by = c("area", "day"))

# S3 method for class 'expjn'
print(x, ...)

# S3 method for class 'expjn'
summary(object, ...)
```

## Arguments

- expo:

  The exposure table: a data frame with the key columns and one or more
  value columns.

- outc:

  The outcome table: a data frame with the key columns.

- by:

  The names of the key columns, in both tables.

- ...:

  Not used.

- object, x:

  An object of class `expjn`.

## Value

An object of class `expjn`, a list: `data`, the joined table (the
outcome table's columns, then the exposure's other columns); `nout`,
`nexp`; `unout`, the keys of the outcome rows left out.
[`summary()`](https://rdrr.io/r/base/summary.html) returns the counts as
a data frame.

## Details

Duplicate keys in either table are refused: a place with two exposures
on one day has no single value to join. Outcome rows without an exposure
are left out of `data` and counted in `nout` (their keys in `unout`);
exposure rows without an outcome are counted in `nexp`. Keys are
compared as text, so a day given as a date in one table must be a date
in the other.

## See also

[`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
[`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md)

## Examples

``` r
ex <- data.frame(
  area = c("a", "a", "b"),
  day = c(1, 2, 1),
  pm = c(10, 12, 8)
)
oc <- data.frame(area = c("a", "a", "b", "b"), day = c(1, 2, 1, 2),
                 crimes = c(3, 1, 0, 2))
j <- expjn(ex, oc)
j
#> expjn: 3 rows joined; 1 outcome row(s) without an exposure left out; 0 exposure row(s) without an outcome
j$data
#>   area day crimes pm
#> 1    a   1      3 10
#> 2    a   2      1 12
#> 3    b   1      0  8
```
