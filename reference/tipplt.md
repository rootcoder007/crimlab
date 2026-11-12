# Tippett plot

For each value on the axis, the share of same-source and of
different-source likelihood ratios above it. A good system keeps the
same-source curve high to the right of 0 and the different-source curve
low; where each crosses 0 shows the rate of misleading evidence.
[`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
records no theorem for this plot.

## Usage

``` r
tipplt(same, diff, ...)
```

## Arguments

- same:

  Likelihood ratios from same-source comparisons: positive numbers.

- diff:

  Likelihood ratios from different-source comparisons: positive numbers.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a list of two data frames, `same` and `diff`, each with
`log10lr` and `share` (the share of ratios at or above it).

## Examples

``` r
tipplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
```
