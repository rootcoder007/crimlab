# Composite index of standardised counts

Ranks units (neighbourhoods, divisions) on several kinds of crime at
once: each category's counts are standardised across the units (minus
the mean, over the standard deviation) and the standardised values are
summed with weights. A category whose counts do not vary contributes 0
to every unit. A unit with no row for a category counts 0 there, and its
standardised value is computed with that 0.

## Usage

``` r
zcomp(x, unit, cat, w = NULL)
```

## Arguments

- x:

  The counts: numbers from 0, one per row.

- unit:

  The unit of each row: labels, at least two units.

- cat:

  The category of each row: labels.

- w:

  `NULL` for equal weights, or the weight of each category: numbers from
  0, named by category or in the order of the categories.

## Value

An object of class `zcomp`, a list: `tab`, a data frame with one row per
unit (`unit`, `comp` the composite, `rank` with 1 the highest and ties
sharing the lowest rank); `cat`, one row per category (`cat`, `mean`,
`sd` and `w`); `z`, the matrix of standardised counts, units by
categories.

## Examples

``` r
d <- data.frame(hood = rep(c("a", "b", "c"), 2),
                kind = rep(c("assault", "theft"), each = 3),
                n = c(5, 9, 1, 30, 22, 41))
zcomp(d$n, d$hood, d$kind)
#> zcomp: composite of 2 standardised categories
#>  unit        comp rank
#>     a -0.10482848    3
#>     b  0.05654365    1
#>     c  0.04828484    2
```
