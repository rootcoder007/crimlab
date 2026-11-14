# Outcome test

The stop rate and the hit rate of a group searched whenever the chance
of guilt exceeds a threshold. Hit rates compare averages, not
thresholds: two groups held to the same threshold can show different hit
rates, and a group held to a lower threshold can show the higher one
(infra-marginality). See `proof("outtst")`.

## Usage

``` r
outtst(q, t, w = NULL)
```

## Arguments

- q:

  The chances of guilt in the group: numbers in the interval 0 to 1.

- t:

  The search threshold: one number.

- w:

  `NULL` for equal weights, or one non-negative weight per element of
  `q` with a positive total; when both are named with unique names, `w`
  is matched to `q` by name.

## Value

A named numeric vector: `stop` (the weighted share searched) and `hit`
(the weighted mean chance of guilt among them; `NaN` when nobody is
searched).

## References

Ayres, I. (2002). Outcome tests of racial disparities in police
practices. *Justice Research and Policy*, 4(1-2), 131-142.
[doi:10.3818/JRP.4.1.2002.131](https://doi.org/10.3818/JRP.4.1.2002.131)

Simoiu, C., Corbett-Davies, S. and Goel, S. (2017). The problem of
infra-marginality in outcome tests for discrimination. *The Annals of
Applied Statistics*, 11(3), 1193-1216.
[doi:10.1214/17-AOAS1058](https://doi.org/10.1214/17-AOAS1058)

## Examples

``` r
outtst(q = c(0.1, 0.5, 0.9), t = 0.3)
#>      stop       hit 
#> 0.6666667 0.7000000 
```
