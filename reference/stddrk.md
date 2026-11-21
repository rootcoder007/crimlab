# Standard drinks in a volume of a beverage

The number of standard drinks in `vol` millilitres of a beverage of
`abv` per cent alcohol by volume: the grams of ethanol (the volume of
ethanol times its density, 0.789 g per mL at 20 degrees Celsius) over
the grams a country calls one drink. Countries differ by a factor of two
and more (Kalinowski and Humphreys 2016), so a count of drinks means
nothing without the standard it was counted in.

## Usage

``` r
stddrk(vol, abv, std = c("ca", "us", "uk", "au", "who"), grams = NULL)
```

## Arguments

- vol:

  The volume, in millilitres.

- abv:

  The alcohol by volume, in per cent (5 for a beer of 5%).

- std:

  The standard: `"ca"` (Canada, 13.45 g, the default), `"us"` (United
  States, 14 g), `"uk"` (a unit, 8 g), `"au"` (Australia, 10 g) or
  `"who"` (World Health Organization, 10 g).

- grams:

  The grams of ethanol in one drink, in place of `std`.

## Value

A numeric vector: standard drinks.

## Details

Errors: `vol` or `abv` not numeric, `vol` negative, `abv` outside 0 to
100; `std` not one of the standards; `grams` not one positive number;
`vol` and `abv` of different lengths (a length of 1 is repeated).

## References

Kalinowski, A. and Humphreys, K. (2016). Governmental standard drink
definitions and low-risk alcohol consumption guidelines in 37 countries.
*Addiction*, 111(7), 1293-1298.
[doi:10.1111/add.13341](https://doi.org/10.1111/add.13341)

## See also

[`ebac()`](https://rootcoder007.github.io/crimlab/reference/ebac.md).

## Examples

``` r
# a 341 mL bottle of 5% beer, a 142 mL glass of 12% wine, 43 mL of 40%
# spirits: one Canadian standard drink each
stddrk(c(341, 142, 43), c(5, 12, 40))
#> [1] 1.0001822 0.9995955 1.0089814
# a 750 mL bottle of 13% wine, by four standards
vapply(c("ca", "us", "uk", "au"), function(s) stddrk(750, 13, s), 0)
#>       ca       us       uk       au 
#> 5.719517 5.494821 9.615938 7.692750 
try(stddrk(341, 500))
#> Error : stddrk: `abv` must be per cent, from 0 to 100.
```
