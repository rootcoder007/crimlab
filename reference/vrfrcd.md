# Prove a recode with a before/after cross-tabulation

Cross-tabulates the original and the recoded column and checks that the
recode is a function of the original label (no category was split
between two new ones), that it follows the declared mapping (a value the
mapping does not name must be unchanged), that no row was lost or added
and that missingness did not change. This is how swapped groups are
caught before any model runs.

## Usage

``` r
vrfrcd(orig, recd, decl)
```

## Arguments

- orig, recd:

  The column before and after the recode.

- decl:

  The declared mapping, a named vector old = new.

## Value

The cross-tabulation as a data frame (`original`, `recoded`, `Freq`),
invisibly; an error at the first failed check.

## See also

[`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
[`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md).

## Examples

``` r
o <- c("M", "F", "F", NA)
m <- c(M = "Male", F = "Female")
vrfrcd(o, grdrcd(o, m), m)
try(vrfrcd(o, c("Female", "Male", "Male", NA), m))
#> Error : vrfrcd: 'F' was mapped to 'Male' but the declared mapping says 'Female'. THIS is how groups get swapped; fix the recode before any model runs.
```
