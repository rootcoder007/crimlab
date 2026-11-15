# Rank correlation and midranks (C backend)

`cspr()` is Spearman's rho, the Pearson correlation of the ranks: it
measures monotone association and is unaffected by any order-preserving
transformation of either variable. `cmidrk()` gives the ranks
themselves; tied values share the average of the ranks they span, which
is what makes the result agree with
[`stats::cor()`](https://rdrr.io/r/stats/cor.html) on tied data.

## Usage

``` r
cspr(x, y)

cmidrk(x)
```

## Arguments

- x, y:

  Numeric (or logical) vectors of the same length.

## Value

`cspr()`: a number in \\\[-1, 1\]\\. `cmidrk()`: the rank of each value,
the length of `x`.

## Details

With \\r_i\\ and \\s_i\\ the midranks of \\x_i\\ and \\y_i\\, \$\$\rho =
\frac{\sum_i (r_i - \bar r)(s_i - \bar s)} {\sqrt{\sum_i (r_i - \bar
r)^2 \sum_i (s_i - \bar s)^2}}\$\$ (Spearman, 1904); the shortcut \\1 -
6 \sum d_i^2 / (n(n^2 - 1))\\ is not used, so ties are handled
correctly. `cspr()` gives `NaN` for an `NA` in either vector and for a
constant vector. `cmidrk()` keeps an `NA` in place and ranks the rest,
as `rank(x, na.last = "keep")` does.

## References

Spearman C (1904). The proof and measurement of association between two
things. *The American Journal of Psychology* 15(1), 72–101.
[doi:10.2307/1412159](https://doi.org/10.2307/1412159)

## See also

[`ccorr()`](https://rootcoder007.github.io/crimlab/reference/cstat.md).

## Examples

``` r
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 9, 16, 25)
cspr(x, y)
#> [1] 1
ccorr(x, y)
#> [1] 0.9737247
xt <- c(1, 2, 2, 2, 5, 5, 7)
yt <- c(3, 1, 1, 4, 4, 9, 2)
all.equal(cspr(xt, yt), stats::cor(xt, yt, method = "spearman"))
#> [1] TRUE
cmidrk(xt)
#> [1] 1.0 3.0 3.0 3.0 5.5 5.5 7.0
cmidrk(c(3, NA, 1, 3))
#> [1] 2.5  NA 1.0 2.5
try(cspr(1:3, 1:4))
#> Error : cspr: `x` and `y` must have the same length.
```
