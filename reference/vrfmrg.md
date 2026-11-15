# Verify recoded counts against the counts a release published

Counts each label of `x` and compares the counts with those the source
published (a table's margins, an annual report). When they differ and
some permutation of the labels would make them agree (up to seven labels
are searched), the labels are attached to the wrong groups, and the
message says which.

## Usage

``` r
vrfmrg(x, pub, tol = 0, strict = TRUE)
```

## Arguments

- x:

  The recoded labels.

- pub:

  The published counts, a named numeric vector (label = count).

- tol:

  The largest difference counted as agreement, one number from 0.

- strict:

  `TRUE` stops at a failure.

## Value

A list of class `vrfmrg` (invisibly when the counts agree): `counts`,
`published`, `ok`, `permutation` (the labels that would fit, or `NULL`)
and `message`.

## Details

A label of `x` that the published counts do not have is a failure. With
`strict` (the default) a failure is an error; otherwise the check is
returned with `ok` `FALSE` and its `message`.

## See also

[`vrfrcd()`](https://rootcoder007.github.io/crimlab/reference/vrfrcd.md),
[`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md).

## Examples

``` r
x <- rep(c("White", "Black", "Indigenous"), c(50, 30, 20))
vrfmrg(x, c(White = 50, Black = 30, Indigenous = 20))
vrfmrg(x, c(White = 30, Black = 50, Indigenous = 20), strict = FALSE)
#> vrfmrg: the counts differ 
#>       label observed published difference
#>       White       50        30         20
#>       Black       30        50        -20
#>  Indigenous       20        20          0
#> vrfmrg: recoded counts do not match the published counts. White:
#> observed 50, published 30; Black: observed 30, published 50;
#> Indigenous: observed 20, published 20. The observed counts match the
#> published ones if the labels are permuted (White -> Black, Black ->
#> White): the labels are attached to the wrong groups.
```
