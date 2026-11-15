# Label the periods of a change table without changing them

A period's value and its name are not the same thing: a fiscal year
keyed by the year it ends (2023 for the year that runs from 2022 into
2023) printed as "2023" names a calendar year the row is not about.
`yoylbl()` adds a column `perlbl` that the writers show in place of the
period; the comparison still runs on the value, so relabelling cannot
change an answer. `fylbl()` names a fiscal year by the years it spans.

## Usage

``` r
yoylbl(x, labels)

fylbl(endyear, sep = "/", short = TRUE)
```

## Arguments

- x:

  A change table
  ([`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)).

- labels:

  A function, a named vector or one label per period.

- endyear:

  The year each fiscal year ends, whole numbers (or text holding them).

- sep:

  The text between the two years.

- short:

  `TRUE` writes the second year with two digits.

## Value

`yoylbl()`: `x` with a column `perlbl`. `fylbl()`: text, `NA` where
`endyear` is missing.

## Details

`labels` is a function of the periods, a vector named by period (a
partial mapping: periods it does not name keep their value), or an
unnamed vector with one label per distinct period, in their order of
appearance.

## See also

[`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md),
[`yoysmry()`](https://rootcoder007.github.io/crimlab/reference/yoysmry.md),
[`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md).

## Examples

``` r
fylbl(c(2023, 2024))
#> [1] "2022/23" "2023/24"
fylbl(2023, sep = "-", short = FALSE)
#> [1] "2022-2023"
# simulated admissions to custody by fiscal year
y <- yoypx(c(120, 131, 98), 2021:2023)
yoylbl(y, fylbl)
#> yoypx: change on the earlier period
#>   prd value prev change     pctchg              flag     pctlo     pctup
#>  2021   120   NA     NA         NA no earlier period        NA        NA
#>  2022   131  120     11   9.166667              <NA> -15.43629 41.033275
#>  2023    98  131    -33 -25.190840              <NA> -43.01728 -2.059398
#>           p       padj signif  perlbl
#>          NA         NA     NA 2020/21
#>  0.52799664 0.52799664  FALSE 2021/22
#>  0.03423381 0.06846762  FALSE 2022/23
yoylbl(y, c("2023" = "the last year"))
#> yoypx: change on the earlier period
#>   prd value prev change     pctchg              flag     pctlo     pctup
#>  2021   120   NA     NA         NA no earlier period        NA        NA
#>  2022   131  120     11   9.166667              <NA> -15.43629 41.033275
#>  2023    98  131    -33 -25.190840              <NA> -43.01728 -2.059398
#>           p       padj signif        perlbl
#>          NA         NA     NA          2021
#>  0.52799664 0.52799664  FALSE          2022
#>  0.03423381 0.06846762  FALSE the last year
```
