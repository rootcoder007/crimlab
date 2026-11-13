# Cumulative risk of a first event by age

The share of a birth cohort that has had a first arrest, conviction or
prison admission by each age, from first events by age and the
population: the hazard at each age is the first events over those still
without one, and the cumulative risk is one minus the product of the
complements, with Greenwood's interval on the log scale of the survival
(Bonczar 2003; Pettit and Western 2004). With death rates the first
events count only among the cohort still alive, and no interval is
given.

## Usage

``` r
cumrsk(age, ev, pop, atrisk = NULL, dth = NULL, level = 0.95)
```

## Arguments

- age:

  The ages, increasing.

- ev:

  The first events at each age: numbers from 0.

- pop:

  The population at each age: positive numbers.

- atrisk:

  `NULL` (the population less those already with an event) or the number
  at risk at each age.

- dth:

  `NULL`, or the death rate at each age, from 0 to below 1.

- level:

  The level of the interval: one number between 0 and 1.

## Value

A data frame of class `cumrsk`: `age`, `hazard`, `atrisk`, `cumrisk`,
`lo`, `hi` (`NA` with death rates).

## References

Pettit, B. and Western, B. (2004). Mass imprisonment and the life
course: race and class inequality in U.S. incarceration. *American
Sociological Review*, 69(2), 151-169.
[doi:10.1177/000312240406900201](https://doi.org/10.1177/000312240406900201)

## Examples

``` r
cumrsk(18:22, c(120, 100, 80, 60, 40), rep(10000, 5))
#> cumrsk: cumulative risk of a first event by age
#>  age      hazard atrisk cumrisk          lo         hi
#>   18 0.012000000  10000   0.012 0.009863582 0.01413181
#>   19 0.010121457   9880   0.022 0.019120830 0.02487072
#>   20 0.008179959   9780   0.030 0.026650783 0.03333769
#>   21 0.006185567   9700   0.036 0.032341858 0.03964431
#>   22 0.004149378   9640   0.040 0.036151577 0.04383306
```
