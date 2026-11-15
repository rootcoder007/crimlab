# Verify a categorical variable that crossed from one program to another

Checks a categorical column after it was exported by one program and
imported by another: the value labels that arrived must agree with the
source's code book, and the decoded values must reproduce the counts the
source reported. When the column carries a `labels` attribute it is
decoded by it
([`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md));
otherwise by `codebook`; otherwise the values are taken as labels.

## Usage

``` r
xfrvrf(imported, soucnts, codebook = NULL, tol = 0, strict = TRUE)
```

## Arguments

- imported:

  The column as imported.

- soucnts:

  The counts the source program reported, a named numeric vector (label
  = count).

- codebook:

  `NULL`, or the source's code book, a named vector code = label.

- tol:

  The largest difference counted as agreement, one number from 0.

- strict:

  `TRUE` stops at a failure.

## Value

A list of class `xfrvrf`: `ok`, `decoded` (a factor), `marginals` (the
[`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md)
check), `codebook` (`TRUE`, `FALSE`, or `NA` when there was nothing to
compare) and `reasons`.

## See also

[`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md),
[`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md).

## Examples

``` r
imp <- structure(c(1, 2, 2, 1, 1), labels = c(No = 1, Yes = 2))
xfrvrf(imp, c(No = 3, Yes = 2), codebook = c("1" = "No", "2" = "Yes"))
#> xfrvrf: as the source reported 
#> code book: agrees 
#> vrfmrg: the counts match 
#>  label observed published difference
#>     No        3         3          0
#>    Yes        2         2          0
xfrvrf(c("No", "Yes", "Yes"), c(No = 1, Yes = 2))
#> xfrvrf: as the source reported 
#> code book: nothing to compare 
#> vrfmrg: the counts match 
#>  label observed published difference
#>     No        1         1          0
#>    Yes        2         2          0
```
