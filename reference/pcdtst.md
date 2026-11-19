# Pesaran's CD test of cross-sectional dependence in a panel

Whether the residuals (or any variable) of a unit-by-period panel move
together across units, the test Andreatta et al. (2023, p. 1376) run on
the residuals of their waste-crime models before deciding a spatial
model is not needed. For each pair of units the correlation `rho_ab` of
their values over the `T_ab` periods they share, about the means over
those periods; `CD = sum sqrt(T_ab) rho_ab /` `sqrt(m)` over the `m`
pairs, which for a balanced panel of `N` units and `T` periods is
`sqrt(2 T / (N (N - 1))) sum rho_ab`, standard normal under no
dependence. Pairs that share fewer than two periods, or along which a
unit does not vary, are left out. This is the statistic of
`plm::pcdtest(test = "cd")`. For dependence between neighbours rather
than across all pairs, use
[`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)
on each period's residuals, as Andreatta et al. also do.

## Usage

``` r
pcdtst(x, unit, time)
```

## Arguments

- x:

  The values of each row: finite numbers (residuals of a panel model,
  typically).

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

## Value

A named numeric vector: `cd`, `p` (two-sided normal), `rho` and `absrho`
(the mean correlation and the mean absolute correlation of the pairs),
`pairs` (`m`), `units` and `periods`.

## References

Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
Digging into waste: an analysis of waste crime in the Italian provinces.
*Regional Studies*, 57(7), 1367-1379.
[doi:10.1080/00343404.2022.2135697](https://doi.org/10.1080/00343404.2022.2135697)

## Examples

``` r
set.seed(3)
d <- expand.grid(year = 2009:2017, prov = paste0("p", 1:20))
common <- rnorm(9)
d$e <- rnorm(180) + 0.8 * common[d$year - 2008]
pcdtst(d$e, d$prov, d$year)
#>           cd            p          rho       absrho        pairs        units 
#> 9.390043e+00 5.997592e-21 2.270751e-01 3.385986e-01 1.900000e+02 2.000000e+01 
#>      periods 
#> 9.000000e+00 
# the residuals of a two-way fit
d$x <- rnorm(180)
d$y <- d$x + d$e
pcdtst(fepnl(d$y, d$x, d$prov, d$year)$resid, d$prov, d$year)
#>           cd            p          rho       absrho        pairs        units 
#>  -2.14401630   0.03203158  -0.05184776   0.29758687 190.00000000  20.00000000 
#>      periods 
#>   9.00000000 
```
