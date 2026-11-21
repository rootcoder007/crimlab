# Estimated blood alcohol concentration (Widmark)

The blood alcohol concentration a person would have after a number of
standard drinks, by Widmark's formula: the alcohol drunk, spread over
the share of the body that takes it up, less what has been eliminated
since drinking began. It is the estimate surveys of drinking report (an
"eBAC") and toxicologists use for back-calculation; it is not a
measurement, and a person's own rate of elimination and body water move
it by a third either way.

## Usage

``` r
ebac(
  drinks,
  weight,
  hours = 0,
  sex = "male",
  height = NULL,
  r = NULL,
  beta = 0.015,
  grams = 14,
  unit = c("kg", "lb")
)
```

## Arguments

- drinks:

  The number of standard drinks.

- weight:

  The body weight, in `unit`.

- hours:

  The hours since drinking began.

- sex:

  `"male"` or `"female"`, one or one per person: which default `r`
  applies. Not used when `r` is given.

- height:

  The height in centimetres, for the factor of Seidl, Jensen and Alt;
  `NULL` (the default) for the fixed factors.

- r:

  The Widmark factor, or `NULL` for the factor of `sex` (and of
  `height`, when given).

- beta:

  The elimination, in grams per 100 mL per hour.

- grams:

  The grams of ethanol in one standard drink: 14 in the United States,
  13.45 in Canada (see
  [`stddrk()`](https://rootcoder007.github.io/crimlab/reference/stddrk.md)),
  13.6 in the survey above.

- unit:

  `"kg"` (the default) or `"lb"`: the unit of `weight`.

## Value

A numeric vector: the concentration in grams per 100 mL, never below 0;
`NA` where an input is missing.

## Details

\$\$BAC = \frac{g \\ d}{10 \\ r \\ W} - \beta t\$\$ in grams per 100 mL,
with `d` drinks of `g` grams of ethanol each, `W` the body weight in
kilograms, `r` the Widmark factor, `beta` the elimination per hour and
`t` the hours since drinking began; a negative value is 0. The form with
fluid ounces and pounds (Matthews and Miller 1979),
`5.14 oz / (r lb) - 0.015 t`, is the same formula: 5.14 is 100 times the
grams in a fluid ounce of ethanol over the grams in a pound.

The default `r`, 0.73 for men and 0.66 for women, are the factors of
that form; Widmark's own means were 0.68 and 0.55. With `height`, `r` is
the factor of Seidl, Jensen and Alt (2000), which follows the body water
of a person of that weight and height:
`0.31608 - 0.004821 W + 0.004632 H` for men and
`0.31223 - 0.006446 W + 0.004466 H` for women, `H` in centimetres. Give
`r` to use another. `beta` 0.015 is the usual mean; reported rates run
from 0.010 to 0.025.

The Canadian Postsecondary Education Alcohol and Drug Use Survey derives
its `ebac_tot` this way (Health Canada 2022, Table 3): the Seidl
factors, 13.6 grams a drink, the drinks and the hours of the heaviest
drinking day, with `beta` 0.017; a value of 0.5 or more is set aside
there as not survivable, and its `ebac_legal` is a value above 0.08. The
guide writes the elimination inside the division by ten,
`(A / (r W) - 0.017 t) / 10`, which is 0.0017 grams per 100 mL an hour;
the values of the public file fall with the hours at that rate, so
`beta = 0.0017` is what reproduces them, and `beta = 0.017` is the rate
of Seidl, Jensen and Alt.

Errors: an argument that is not numeric; a negative `drinks` or `hours`;
a `weight`, `height` or `r` that is not positive (a Seidl factor that is
not positive too: a weight and height no person has); `beta` not one
number from 0; `grams` not one positive number; a `sex` other than
`"male"` or `"female"`; arguments of different lengths (a length of 1 is
repeated).

## References

Widmark, E. M. P. (1932). *Die theoretischen Grundlagen und die
praktische Verwendbarkeit der gerichtlich-medizinischen
Alkoholbestimmung*. Urban & Schwarzenberg, Berlin.

Matthews, D. B. and Miller, W. R. (1979). Estimating blood alcohol
concentration: two computer programs and their applications in therapy
and research. *Addictive Behaviors*, 4(1), 55-60.
[doi:10.1016/0306-4603(79)90021-2](https://doi.org/10.1016/0306-4603%2879%2990021-2)

Seidl, S., Jensen, U. and Alt, A. (2000). The calculation of blood
ethanol concentrations in males and females. *International Journal of
Legal Medicine*, 114(1-2), 71-77.
[doi:10.1007/s004140000154](https://doi.org/10.1007/s004140000154)

Health Canada (2022). *Canadian Postsecondary Education Alcohol and Drug
Use Survey, 2021-22: Public Use Microdata File user guide*. Table 3,
derived variables `ebac_tot` and `ebac_legal`.

## See also

[`stddrk()`](https://rootcoder007.github.io/crimlab/reference/stddrk.md)
for standard drinks,
[`hvydrk()`](https://rootcoder007.github.io/crimlab/reference/hvydrk.md)
and
[`lrdg()`](https://rootcoder007.github.io/crimlab/reference/hvydrk.md)
for heavy drinking and the low-risk guideline;
[`etohlr()`](https://rootcoder007.github.io/crimlab/reference/etohlr.md)
for ethanol found after death.

## Examples

``` r
# four drinks in two hours, 180 lb
ebac(4, 180, 2, unit = "lb")
#> [1] 0.06395652
# the same person at 0, 1, 2 and 3 hours, and a woman of the same weight
ebac(4, 81.6, 0:3)
#> [1] 0.09401021 0.07901021 0.06401021 0.04901021
ebac(4, 81.6, 2, sex = "female")
#> [1] 0.07398099
# the survey's way: body water from weight and height, 13.6 g a drink
ebac(5, 62, 3, sex = "female", height = 164, beta = 0.017, grams = 13.6)
#> [1] 0.119042
# over 0.08?
ebac(c(2, 5), c(70, 70), 1) >= 0.08
#> [1] FALSE  TRUE
try(ebac(4, 0, 2))
#> Error : ebac: `weight` must be positive numbers.
```
