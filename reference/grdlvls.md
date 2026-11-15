# A factor with declared, verified levels

Builds a factor whose levels are exactly the ones declared, in the
declared order, and refuses any value outside them. A declared reference
must be the first level, so the group every model contrast is measured
against is stated rather than left to alphabetical order.

## Usage

``` r
grdlvls(x, lvls, ref = NULL)
```

## Arguments

- x:

  The labels: a vector or a factor.

- lvls:

  The levels, a non-empty character vector without `NA`, the reference
  first.

- ref:

  `NULL`, or the reference level, which must be `lvls[1]`.

## Value

A factor with levels `lvls`.

## See also

[`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md).

## Examples

``` r
grdlvls(c("Male", "Female", "Male"), c("Male", "Female"), ref = "Male")
#> [1] Male   Female Male  
#> Levels: Male Female
try(grdlvls(c("Male", "Female", "X"), c("Male", "Female")))
#> Error : grdlvls: values outside the declared levels: 'X'.
try(grdlvls("Male", c("Male", "Female"), ref = "Female"))
#> Error : grdlvls: the declared reference 'Female' is not the first level ('Male'); put the reference first in `lvls`.
```
