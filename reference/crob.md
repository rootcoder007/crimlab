# Quantiles, median and robust spread (C backend)

`cqnt()` is the type-7 quantile, R's default, so it agrees with
`stats::quantile(x, probs)`. `cmed()` is the 50% point. `cmad()` is the
median absolute deviation scaled by `const` so that it estimates the
standard deviation of a normal sample. `ciqr()` is the interquartile
range and `ctkyfnc()` the outlier fences drawn `k` IQRs beyond the
quartiles. A single corrupted row can move a mean or a variance
arbitrarily far, but moves a median or a MAD hardly at all, which is
what is wanted when screening a column of a released table for values to
inspect by hand.

## Usage

``` r
cqnt(x, probs = c(0, 0.25, 0.5, 0.75, 1))

cmed(x)

cmad(x, const = 1.4826)

ciqr(x)

ctkyfnc(x, k = 1.5)
```

## Arguments

- x:

  Numeric (or logical) vector.

- probs:

  Probabilities, each in \\\[0, 1\]\\.

- const:

  The scale factor of `cmad()`, one number from 0: the default 1.4826 is
  the rounded constant
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html) uses, so the two
  agree exactly; `1` gives the unscaled median deviation.

- k:

  The fence width in IQRs, one number from 0 (1.5 is Tukey's convention,
  3 the usual far-out cut).

## Value

`cqnt()`: a numeric vector the length of `probs`, named as
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) names it.
`cmed()`, `cmad()`, `ciqr()`: a number. `ctkyfnc()`: a named numeric
vector, `lower` and `upper`.

## Details

The type-7 quantile at probability \\p\\ interpolates between the order
statistics \\x\_{(j)}\\ and \\x\_{(j+1)}\\ with \\h = (n - 1)p + 1\\,
\\j = \lfloor h \rfloor\\: \\Q(p) = x\_{(j)} + (h - j) (x\_{(j+1)} -
x\_{(j)})\\ (Hyndman and Fan, 1996). The MAD is \\c \cdot
\mathrm{median}\_i \|x_i - \mathrm{median}(x)\|\\; the IQR \\Q(0.75) -
Q(0.25)\\; the fences \\Q(0.25) - k \cdot \mathrm{IQR}\\ and \\Q(0.75) +
k \cdot \mathrm{IQR}\\ (Tukey, 1977).

Missing values are not removed: an `NA` or `NaN` anywhere makes every
result `NaN` (where
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) stops). An
empty `x` gives `NaN`. `cqnt()` refuses an empty `probs` and any
probability outside \\\[0, 1\]\\ or missing.

## References

Hyndman RJ, Fan Y (1996). Sample quantiles in statistical packages. *The
American Statistician* 50(4), 361–365.
[doi:10.1080/00031305.1996.10473566](https://doi.org/10.1080/00031305.1996.10473566)

Tukey JW (1977). *Exploratory Data Analysis*. Addison-Wesley.

## See also

[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md),
[`ctrmmean()`](https://rootcoder007.github.io/crimlab/reference/ctrm.md).

## Examples

``` r
x <- c(2, 4, 4, 4, 5, 5, 7, 9)
cqnt(x, c(0.25, 0.5, 0.75))
#> 25% 50% 75% 
#> 4.0 4.5 5.5 
cmed(x)
#> [1] 4.5
cmad(x)
#> [1] 0.7413
all.equal(cmad(x), stats::mad(x))
#> [1] TRUE
ciqr(x)
#> [1] 1.5
ctkyfnc(x)
#> lower upper 
#>  1.75  7.75 
wild <- c(x, 1000)
f <- ctkyfnc(wild)
wild[wild < f[["lower"]] | wild > f[["upper"]]]
#> [1] 1000
cmed(c(1, NA, 3))
#> [1] NaN
try(cqnt(x, 1.5))
#> Error : cqnt: `probs` must lie in [0, 1].
```
