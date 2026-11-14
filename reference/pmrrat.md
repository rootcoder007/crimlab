# Postmortem redistribution ratio

After death drugs move from organs into the central blood vessels, so
that heart blood can read well above the concentration at death. The
ratio of the central (heart) to the peripheral (femoral) concentration
flags it: at most `cut[1]` (1) minimal, at most `cut[2]` (2) modest,
above that significant, when the peripheral value is the better guide to
the concentration at death.

## Usage

``` r
pmrrat(cen, per, cut = c(1, 2))
```

## Arguments

- cen, per:

  The central and peripheral blood concentrations: positive numbers, one
  per case.

- cut:

  The two cut points of the ratio: positive, increasing.

## Value

A data frame with one row per case: `cen`, `per`, `ratio` (`cen / per`)
and `redist`, a factor (`minimal`, `modest`, `significant`).

## References

Yarema, M. C. and Becker, C. E. (2005). Key concepts in postmortem drug
redistribution. *Clinical Toxicology*, 43(4), 235-241.
[doi:10.1081/CLT-58950](https://doi.org/10.1081/CLT-58950)

## Examples

``` r
pmrrat(c(0.4, 1.5, 0.9), c(0.5, 1, 0.2))
#>   cen per ratio      redist
#> 1 0.4 0.5   0.8     minimal
#> 2 1.5 1.0   1.5      modest
#> 3 0.9 0.2   4.5 significant
```
