# Empirical cross-entropy plot

The information a set of likelihood ratios still leaves missing, across
prior odds: for the ratios as reported, for the same ratios after
calibration, and for a system that always reports 1. Reported ratios
above the last curve do harm at that prior. At even prior odds the first
two curves pass through the values of
[`cllrmin()`](https://rootcoder007.github.io/crimlab/reference/cllrmin.md).
[`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
records no theorem for this plot.

## Usage

``` r
eceplt(same, diff, lo = seq(-2.5, 2.5, by = 0.05), ...)
```

## Arguments

- same:

  Likelihood ratios from same-source comparisons: positive numbers.

- diff:

  Likelihood ratios from different-source comparisons: positive numbers.

- lo:

  The prior log10 odds to draw at.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a data frame: `lo`, `ece`, `min` and `ref`.

## References

Ramos, D. and Gonzalez-Rodriguez, J. (2013). Reliable support: measuring
calibration of likelihood ratios. *Forensic Science International*,
230(1-3), 156-169.
[doi:10.1016/j.forsciint.2013.04.014](https://doi.org/10.1016/j.forsciint.2013.04.014)

## Examples

``` r
eceplt(same = c(50, 200, 8, 0.5), diff = c(0.01, 0.2, 1.5, 0.001))
```
