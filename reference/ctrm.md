# Trimmed and winsorized means (C backend)

Two ways to stop a handful of extreme rows dominating a column's centre.
`ctrmmean()` DISCARDS the `floor(n * trim)` largest and smallest values,
as `mean(x, trim = )` does. `cwinsmean()` instead PULLS THEM IN to the
most extreme surviving values, so every observation still counts,
usually the better choice when the extremes are real measurements rather
than errors.

## Usage

``` r
ctrmmean(x, trim = 0.1)

cwinsmean(x, trim = 0.1)
```

## Arguments

- x:

  Numeric (or logical) vector.

- trim:

  The proportion trimmed (or winsorized) at each end, one number from 0
  to 0.5.

## Value

A number in the units of `x`.

## Details

With \\g = \lfloor n \cdot \mathrm{trim} \rfloor\\ and the sorted values
\\x\_{(1)} \le \dots \le x\_{(n)}\\, the trimmed mean is \\\frac{1}{n -
2g} \sum\_{i=g+1}^{n-g} x\_{(i)}\\ and the winsorized mean replaces
\\x\_{(1)}, \dots, x\_{(g)}\\ by \\x\_{(g+1)}\\ and \\x\_{(n-g+1)},
\dots, x\_{(n)}\\ by \\x\_{(n-g)}\\ before averaging all \\n\\ values.
At `trim = 0.5` both return the median, as
[`mean()`](https://rdrr.io/r/base/mean.html) does. An empty `x`, or any
`NA` or `NaN`, gives `NaN`.

## See also

[`cmed()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
[`cmad()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
[`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md).

## Examples

``` r
x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
mean(x)
#> [1] 14.5
ctrmmean(x, 0.1)
#> [1] 5.5
cwinsmean(x, 0.1)
#> [1] 5.5
all.equal(ctrmmean(x, 0.2), mean(x, trim = 0.2))
#> [1] TRUE
all.equal(ctrmmean(x, 0.5), cmed(x))
#> [1] TRUE
try(cwinsmean(x, 0.7))
#> Error : cwinsmean: `trim` must be one number from 0 to 0.5.
```
