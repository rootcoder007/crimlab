# Risks behind an odds ratio

An odds ratio is not a relative risk: it lies further from 1, and by
more the more common the outcome. Given the odds ratio, the overall rate
of the outcome and the share exposed, this returns the two risks and the
relative risk they imply. See `proof("orrisk")`.

## Usage

``` r
orrisk(or, rate, share)
```

## Arguments

- or:

  The odds ratio: one positive number.

- rate:

  The overall rate of the outcome, strictly between 0 and 1.

- share:

  The exposed share of the population, strictly between 0 and 1.

## Value

A named numeric vector: `exposed` and `unexposed` (the two risks), `rr`
(their ratio) and `or`. A risk closer to 1 than a double can hold is
returned as 1.

## Examples

``` r
orrisk(or = 3, rate = 0.4, share = 0.5)
#>   exposed unexposed        rr        or 
#> 0.5282202 0.2717798 1.9435596 3.0000000 
```
