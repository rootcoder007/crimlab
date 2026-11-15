# Correlation of counts per unit between categories

Counts the events of each category (offence type, call type) by unit
(neighbourhood, beat, district), over the union of the units with a 0
where a category has none, and gives the Pearson correlations of those
counts between every pair of categories: whether the same places score
high on both.

## Usage

``` r
catcor(units)
```

## Arguments

- units:

  A named list, one element per category: the unit of each event (`NA`
  dropped).

## Value

A list of class `catcor`: `r` (categories by categories, `NA` for a
category of one count everywhere), `counts` (units by categories),
`units` and `cats` (sorted).

## References

Pearson, K. (1896). Mathematical contributions to the theory of
evolution. III. Regression, heredity, and panmixia. *Philosophical
Transactions of the Royal Society A*, 187, 253-318.

## See also

Other counts over time:
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
[`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md),
[`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md),
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
[`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)

## Examples

``` r
# simulated events of three categories in 20 areas
set.seed(1)
a <- sprintf("%02d", 1:20)
u <- list(theft = sample(a, 200, TRUE), assault = sample(a, 80, TRUE),
          fraud = sample(a, 50, TRUE))
catcor(u)$r
#>             assault       fraud       theft
#> assault  1.00000000 -0.39207576  0.06201737
#> fraud   -0.39207576  1.00000000 -0.05721296
#> theft    0.06201737 -0.05721296  1.00000000
```
