# Summary statistics accumulated in blocks

An accumulator for the count, mean, variance, skewness and kurtosis of a
numeric column that does not fit in memory, or that arrives in pieces (a
file read in chunks, results from several machines). Start one with
`onlsmry()`, feed it blocks with `smryupd()`, combine two with
`smrymrg()` and read the result with `smrystat()`.

## Usage

``` r
onlsmry(x = NULL)

smryupd(acc, x)

smrymrg(a, b)

smrystat(acc)

# S3 method for class 'onlsmry'
print(x, ...)
```

## Arguments

- x:

  Numeric vector (finite; `NA` is skipped). Text that reads as numbers
  is taken as those numbers.

- acc, a, b:

  Accumulators from `onlsmry()`.

- ...:

  Not used.

## Value

`onlsmry()`, `smryupd()` and `smrymrg()`: an accumulator, a list of
class `onlsmry` with `state`. `smrystat()`: a named numeric vector `n`,
`mean`, `variance`, `sd`, `skewness`, `kurtosis`.

## Details

The merge is exact algebra, not an approximation: it uses Chan, Golub
and LeVeque's parallel combination of central sums, extended to the
third and fourth moments by Terriberry, so accumulating a column in
blocks gives the answer of one pass over the whole column, to the
rounding of double arithmetic. Each block's sums are taken in two
passes, as [`var()`](https://rdrr.io/r/stats/cor.html) takes them, so
values far from zero keep their digits. Memory stays constant in the
number of blocks.

The state is \\(n, \bar{x}, M_2, M_3, M_4)\\ with \\M_k = \sum (x_i -
\bar{x})^k\\. `smrystat()` reports the sample variance \\M_2 / (n -
1)\\, the moment skewness \\g_1 = (M_3/n) / (M_2/n)^{3/2}\\ and the
excess kurtosis \\g_2 = (M_4/n) / (M_2/n)^2 - 3\\ (the biased moment
estimators, without small-sample correction). `NA` values are dropped
from every block. Undefined quantities are `NaN`: the mean for \\n =
0\\, the variance for \\n \< 2\\, the skewness for \\n \< 3\\ and the
kurtosis for \\n \< 4\\, and both shape statistics for a constant
column.

The object is a plain list and is copied on assignment like any other R
value, so `smryupd()` returns the new accumulator, which must be kept:
`acc <- smryupd(acc, block)`.

Errors: a block that is not numeric or holds an infinite value; an
accumulator that does not come from `onlsmry()`.

## References

Chan, T. F., Golub, G. H. and LeVeque, R. J. (1983). Algorithms for
computing the sample variance: analysis and recommendations. The
American Statistician 37(3), 242-247.
[doi:10.1080/00031305.1983.10483115](https://doi.org/10.1080/00031305.1983.10483115)

## Examples

``` r
set.seed(1)
x <- stats::rnorm(1000) # simulated
# accumulate in blocks of 100
acc <- onlsmry()
for (i in seq(1, 1000, by = 100)) acc <- smryupd(acc, x[i:(i + 99)])
smrystat(acc)
#>             n          mean      variance            sd      skewness 
#>  1.000000e+03 -1.164814e-02  1.071051e+00  1.034916e+00 -1.916710e-02 
#>      kurtosis 
#> -1.775465e-03 
all.equal(smrystat(acc)[["variance"]], stats::var(x))
#> [1] TRUE
# two accumulators built independently can be merged
a <- smryupd(onlsmry(), x[1:400])
b <- smryupd(onlsmry(), x[401:1000])
all.equal(smrystat(smrymrg(a, b)), smrystat(acc))
#> [1] TRUE
# an empty accumulator reports nothing rather than zero
smrystat(onlsmry())
#>        n     mean variance       sd skewness kurtosis 
#>        0      NaN      NaN      NaN      NaN      NaN 
try(onlsmry(c("a", "b")))
#> Error : onlsmry: `x` must be numeric, not text such as "a".
```
