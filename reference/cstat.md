# Mean, variance and correlation (C backend)

`cmn()`, `cvari()` and `ccorr()` compute the mean, the sample variance
and the Pearson correlation in compiled code, the same kernels that
crimlab's other front ends (the command line, Python and the bindings)
call, so a result computed here is the result they give. For everyday
work [`mean()`](https://rdrr.io/r/base/mean.html),
[`var()`](https://rdrr.io/r/stats/cor.html) and
[`cor()`](https://rdrr.io/r/stats/cor.html) give the same numbers.

## Usage

``` r
cmn(x)

cvari(x)

ccorr(x, y)
```

## Arguments

- x, y:

  Numeric (or logical) vectors; `x` and `y` of the same length in
  `ccorr()`.

## Value

A number.

## Details

With \\n\\ values, \$\$\bar x = \frac{1}{n}\sum_i x_i, \qquad s^2 =
\frac{1}{n - 1}\sum_i (x_i - \bar x)^2, \qquad r = \frac{\sum_i (x_i -
\bar x)(y_i - \bar y)} {\sqrt{\sum_i (x_i - \bar x)^2 \sum_i (y_i - \bar
y)^2}}.\$\$ The variance uses the \\n - 1\\ (sample) denominator, as
[`stats::var()`](https://rdrr.io/r/stats/cor.html) does. The mean is
base R's (an extended-precision sum and one corrective pass), and a sum
of finite values that would overflow a double (two values of `1e308`)
still gives the finite mean. The correlation is computed from centred
sums and clamped to \\\[-1, 1\]\\.

There is no `na.rm`: an `NA` or `NaN` gives `NA` or `NaN`; remove
missing values first, for example with
[`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html). Edge cases
return `NaN` rather than an error: an empty vector for the mean, fewer
than two values for the variance, a constant vector for the correlation.
Text, a factor, a list or a data frame is refused; logical input is read
as 0 and 1. `ccorr()` refuses `x` and `y` of different lengths.

## See also

[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md),
[`cspr()`](https://rootcoder007.github.io/crimlab/reference/crank.md),
[`ccovm()`](https://rootcoder007.github.io/crimlab/reference/ccovm.md).

## Examples

``` r
cmn(1:10)
#> [1] 5.5
cmn(c(1, 2, NA))
#> [1] NA
cmn(c(1e308, 1e308))
#> [1] 1e+308
cvari(c(2, 4, 4, 4, 5, 5, 7, 9))
#> [1] 4.571429
all.equal(cvari(1:10), stats::var(1:10))
#> [1] TRUE
ccorr(1:10, (1:10)^2)
#> [1] 0.9745586
ccorr(1:10, 10:1)
#> [1] -1
try(cmn(c("1", "2")))
#> Error : cmn: `x` must be numbers (or TRUE and FALSE).
try(ccorr(1:3, 1:4))
#> Error : ccorr: `x` and `y` must have the same length.
```
