# Refuse a categorical column where a 0/1 indicator is required

Checks that a column meant as a numeric 0/1 indicator (a treatment, an
outcome) is one. A factor or text column is refused rather than
converted: [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a
factor returns the level indices (1, 2, ...), not the data, and a
mis-ordered level relabels every observation. Encode it explicitly
first, for example `as.integer(x == "treated")`.

## Usage

``` r
grdbin(x, col)
```

## Arguments

- x:

  The column.

- col:

  Its name, for the message.

## Value

`TRUE`, invisibly; an error otherwise.

## See also

[`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md).

## Examples

``` r
grdbin(c(0, 1, 1, NA), "custody")
try(grdbin(factor(c("no", "yes")), "custody"))
#> Error : grdbin: column 'custody' is categorical ('no', 'yes'...). It is not coerced: as.numeric() on a factor returns level indices, not the data. Encode it first, for example as.integer(x == "treated").
try(grdbin(c(0, 1, 2), "custody"))
#> Error : grdbin: column 'custody' must be 0 or 1 (saw: 0, 1, 2).
```
