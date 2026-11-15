# Lotka-Volterra rates of a yearly count series

Treats yearly counts as the prey \\x(t)\\ and their three-year trailing
mean as a stand-in predator \\y(t)\\ (an enforcement series can take its
place) in \\\dot x = \alpha x - \beta x y\\, \\\dot y = \delta x y -
\gamma y\\: \\\alpha\\ is the median log growth of \\x\\ plus 0.05,
\\\gamma\\ minus the median log growth of \\y\\ plus 0.05, \\\beta =
\alpha / \mathrm{median}(y)\\, \\\delta = \gamma / \mathrm{median}(x)\\,
and the linearised cycle has period \\2\pi / \sqrt{\alpha\gamma}\\
(D'Orsogna and Perc 2015, section 3.4). The rates describe the series;
they are not a fitted dynamical model.

## Usage

``` r
lvfit(y, year = NULL)
```

## Arguments

- y:

  The yearly counts (counts below 1 count as 1), five at least.

- year:

  `NULL` or the year of each count.

## Value

A list of class `lvfit`: `alpha`, `beta`, `gamma`, `delta`, `period`,
`year`, `prey` and `pred`.

## References

D'Orsogna, M. R. and Perc, M. (2015). Statistical physics of crime: a
review. *Physics of Life Reviews*, 12, 1-21.
[doi:10.1016/j.plrev.2014.11.001](https://doi.org/10.1016/j.plrev.2014.11.001)

## See also

Other counts over time:
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
[`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md),
[`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md),
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
[`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)

## Examples

``` r
lvfit(c(410, 455, 430, 470, 520, 498, 515), year = 2017:2023)$period
#> [1] 299.0938
```
