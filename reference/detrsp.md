# Offending chosen under certainty and severity of sanction

An offender picks the level of offending `x` from a menu to maximise the
benefit less the expected sanction, `ben - p * sanc`, with the sanction
increasing in `x`. For any shape of the benefit (no concavity, no
smoothness) the chosen level never rises with the certainty `p`, and
never rises when the sanction schedule grows faster: convexity buys a
unique and smooth response, not its direction. Gives the chosen level at
each certainty; of tied optima (within `1e-12`), the largest.

## Usage

``` r
detrsp(x, ben, sanc, p)
```

## Arguments

- x:

  The levels of offending on the menu: finite numbers, no level twice
  (any order).

- ben:

  The benefit of each level: finite numbers.

- sanc:

  The sanction at each level: finite numbers, increasing with `x`.

- p:

  The certainties of sanction to evaluate: numbers from 0.

## Value

A data frame with one row per certainty: `p`, `xopt` (the level chosen)
and `value` (its net benefit).

## References

Becker, G. S. (1968). Crime and punishment: an economic approach.
*Journal of Political Economy*, 76(2), 169-217.
[doi:10.1086/259394](https://doi.org/10.1086/259394)

## Examples

``` r
x <- 0:10
detrsp(x, ben = 3 * sqrt(x) - (x %% 3 == 0), sanc = x^1.5,
       p = c(0.1, 0.3, 0.5, 1))
#>     p xopt    value
#> 1 0.1   10 6.324555
#> 2 0.3    4 3.600000
#> 3 0.5    2 2.828427
#> 4 1.0    1 2.000000
```
