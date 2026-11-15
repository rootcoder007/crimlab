# Decode a labelled import by its value labels

Decodes a column of codes by the value labels it carries (the `labels`
attribute that readers of SPSS, Stata and SAS files attach, a named
vector label = code), or by `vlabs` (code = label) when it carries none.
The codes alone are not the categories: decoding by code keeps every
observation with its own label whatever order the labels are listed in.

## Usage

``` r
declbld(x, vlabs = NULL)
```

## Arguments

- x:

  The codes, with or without a `labels` attribute.

- vlabs:

  `NULL`, or a named character vector code = label, used when `x` has no
  `labels` attribute.

## Value

A factor.

## Details

The levels of the result are the labels in the order of their codes
(numerically when the codes are numbers). A column with neither a
`labels` attribute nor `vlabs` is refused.

## See also

[`deccodes()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
[`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md).

## Examples

``` r
x <- structure(c(2, 1, 3, 2), labels = c(Female = 1, Male = 2, Other = 3))
declbld(x)
#> [1] Male   Female Other  Male  
#> attr(,"recode")
#> attr(,"recode")$mapping
#>      1      2      3 
#> Female   Male  Other 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
#> Levels: Female Male Other
declbld(c(10, 20, 10), c("20" = "Remand", "10" = "Sentenced"))
#> [1] Sentenced Remand    Sentenced
#> attr(,"recode")
#> attr(,"recode")$mapping
#>        20        10 
#>    Remand Sentenced 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
#> Levels: Sentenced Remand
try(declbld(c(1, 2)))
#> Error : declbld: `x` carries no `labels` attribute and `vlabs` was not given; the codes alone are not the categories.
```
