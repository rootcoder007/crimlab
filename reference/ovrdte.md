# Dates written out in reports, as ISO dates

Finds the first date written in each text, in English or French, and
gives it as `YYYY-MM-DD`, the form that sorts and that
[`as.Date()`](https://rdrr.io/r/base/as.Date.html) reads.
`"January 5, 2023"`, `"Jan. 5th 2023"`, `"5 January 2023"`,
`"5th of January, 2023"`, `"5 janvier 2023"` and `"2023-01-05"` all give
`"2023-01-05"`.

## Usage

``` r
ovrdte(x)
```

## Arguments

- x:

  The texts: a character vector.

## Value

A character vector the length of `x`: the first date of each as
`YYYY-MM-DD`, `NA` where there is none (or `x` is `NA`).

## Details

Month names are English or French, whole or abbreviated (with or without
a full stop), in any case, French accents optional; the day may carry an
ordinal (`"1st"`, `"1er"`). A day that does not exist
(`"February 30, 2023"`) is passed over, and the search goes on to the
next date in the text. The reading never depends on the session's
locale.

## Examples

``` r
ovrdte(c("January 5, 2023", "le 1er février 2021", "no date"))
#> [1] "2023-01-05" "2021-02-01" NA          
ovrdte("On February 30, 2023 or rather March 2, 2023, the SIU ...")
#> [1] "2023-03-02"
as.Date(ovrdte("Date: April 28, 2023"))
#> [1] "2023-04-28"
```
