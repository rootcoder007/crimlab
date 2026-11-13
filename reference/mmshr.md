# Sentences at, below and above a mandatory minimum

The shares of sentences at the mandatory minimum, below it (where a
court found a way around it) and above it, and the mean excess of those
above, overall and by group (offence, court).

## Usage

``` r
mmshr(sent, min, grp = NULL)
```

## Arguments

- sent:

  The sentence of each case: numbers.

- min:

  The mandatory minimum: one number, or one per case.

- grp:

  `NULL`, or the group of each case: labels.

## Value

A data frame of class `mmshr` with one row per group and a last row
`"all"`: `group`, `n`, `at`, `below`, `above` (shares) and `mexcess` (0
when none is above).

## Examples

``` r
mmshr(c(365, 400, 300, 365, 730), 365, c("a", "a", "b", "b", "b"))
#> mmshr: sentences against the minimum
#>  group n        at     below     above mexcess
#>      a 2 0.5000000 0.0000000 0.5000000      35
#>      b 3 0.3333333 0.3333333 0.3333333     365
#>    all 5 0.4000000 0.2000000 0.4000000     200
```
