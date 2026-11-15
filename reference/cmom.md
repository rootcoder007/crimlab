# Mean, variance, skewness and kurtosis (C backend)

The four standard descriptive moments of a numeric vector in one call:
where it is centred, how spread out it is, whether it leans to one side
and how heavy its tails are.

## Usage

``` r
cmom(x)
```

## Arguments

- x:

  Numeric (or logical) vector.

## Value

A named numeric vector: `mean`, `variance`, `skewness` and `kurtosis`.

## Details

The mean is base R's (an extended-precision sum and a corrective pass);
the deviations from it are divided by their largest magnitude before any
power is taken, and the scale is restored at the end, which avoids both
the cancellation of the textbook formula and the overflow of fourth
powers beyond about \\10^{77}\\. The variance uses the `n - 1`
denominator of [`stats::var()`](https://rdrr.io/r/stats/cor.html). The
shape statistics are the sample moment ratios \$\$g_1 = m_3 / m_2^{3/2},
\qquad g_2 = m_4 / m_2^2 - 3,\$\$ with \\m_k = n^{-1} \sum_i (x_i - \bar
x)^k\\, so kurtosis is EXCESS kurtosis (a normal sample sits near zero).
These are the uncorrected ratios (type 1 of Joanes and Gill, 1998).

Edge cases are `NaN`, never an error: an empty vector gives four `NaN`;
one value gives its mean and `NaN` for the rest; skewness needs three
values and kurtosis four; any `NA` or `NaN` makes every entry `NaN`. A
constant vector has variance 0 and `NaN` shape statistics.

## References

Joanes DN, Gill CA (1998). Comparing measures of sample skewness and
kurtosis. *Journal of the Royal Statistical Society: Series D (The
Statistician)* 47(1), 183–189.
[doi:10.1111/1467-9884.00122](https://doi.org/10.1111/1467-9884.00122)

## See also

[`csdev()`](https://rootcoder007.github.io/crimlab/reference/cspread.md),
[`cqnt()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
[`cmad()`](https://rootcoder007.github.io/crimlab/reference/crob.md).

## Examples

``` r
cmom(c(2, 4, 4, 4, 5, 5, 7, 9))
#>      mean  variance  skewness  kurtosis 
#>  5.000000  4.571429  0.656250 -0.218750 
m <- cmom(1:10)
all.equal(m[["variance"]], stats::var(1:10))
#> [1] TRUE
cmom(c(-2, -1, 0, 1, 2))[["skewness"]]
#> [1] 0
cmom(c(1, 2))
#>     mean variance skewness kurtosis 
#>      1.5      0.5      NaN      NaN 
```
