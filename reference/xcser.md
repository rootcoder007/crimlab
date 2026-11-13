# Correlation of two series, at lags

Correlates two series period by period: Pearson's r with Fisher's
interval and the t-test p-value of
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html), and the
cross-correlations at lags \\-K..K\\ as
[`stats::ccf()`](https://rdrr.io/r/stats/acf.html) estimates them (lag
`k` is the correlation of `x[t + k]` with `y[t]`; the default `K` is
[`ccf()`](https://rdrr.io/r/stats/acf.html)'s,
`floor(10 log10(n / 2))`). Two named series are matched by name (the
periods both have, in the order of `x`), so counts from two sources for
the same years line up whatever years each covers. With `dif` the
changes from one period to the next are correlated instead of the
levels: two series that both trend correlate in level whether or not
they move together. Every correlation lies in \\\[-1, 1\]\\
(Cauchy-Schwarz, `proof("xcser")`).

## Usage

``` r
xcser(x, y, lagmax = NULL, dif = FALSE, level = 0.95)
```

## Arguments

- x, y:

  The series: finite numbers, of one length, or both named by period.

- lagmax:

  `NULL` or the largest lag: one whole number from 0, below the number
  of periods.

- dif:

  `TRUE` to correlate the period-to-period changes.

- level:

  The confidence level: one number strictly between 0 and 1.

## Value

An object of class `xcser`, a list: `cor` (`r`, `lo`, `hi`, `p`, `n`);
`ccf` (`lag`, `r`); `common`, the periods matched by name (`NULL`
without names); and `level`.

## Examples

``` r
a <- c("2019" = 120, "2020" = 98, "2021" = 105, "2022" = 130, "2023" = 141)
b <- c("2018" = 50, "2019" = 61, "2020" = 44, "2021" = 52, "2022" = 66,
       "2023" = 70)
xcser(a, b)
#> xcser: 5 periods; r 0.9840288 (95% 0.7719537 to 0.9989935), p 0.002417117
#>  lag          r
#>   -3 -0.3061542
#>   -2 -0.4281115
#>   -1  0.1949278
#>    0  0.9840288
#>    1  0.3245041
#>    2 -0.4603044
#>    3 -0.3987087
xcser(a, b, dif = TRUE)$cor
#>           r        lo        hi          p n
#> 1 0.9734683 0.1922192 0.9994666 0.02653166 4
```
