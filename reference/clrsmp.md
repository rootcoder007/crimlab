# Clearance-rate decomposition and Simpson's paradox

Splits the change in an aggregate clearance rate between two periods,
areas or services (A and B) into a rate effect, from each offence type
clearing more or less, and a mix effect, from the workload moving
between types that clear at different rates (Kitagawa's decomposition).
The two always add up to the change. When every type clears better in B
and the aggregate still falls, the fall is all mix: Simpson's paradox,
and the reason an aggregate clearance rate cannot be read as performance
without the split by type. See `proof("clrsmp")`.

## Usage

``` r
clrsmp(clra, reca, clrb, recb)
```

## Arguments

- clra:

  Offences cleared in A, by type: non-negative numbers, as a numeric
  vector, a 1-d table, a list of single numbers or a one-row data frame.

- reca:

  Offences recorded in A, by type: positive numbers, in the forms of
  `clra`; named inputs are matched to `clra` by name, unnamed ones by
  position (only names protect against a mismatch).

- clrb:

  Offences cleared in B, by type: non-negative numbers, in the forms of
  `clra`, matched by name.

- recb:

  Offences recorded in B, by type: positive numbers, in the forms of
  `clra`, matched by name.

## Value

A list. `parts`, a data frame with one row per type, in the order of
`clra`: `part`, `ratea`, `rateb` (clearance rates), `sharea`, `shareb`
(shares of the recorded workload), `rateeff` and `mixeff`. `total`, a
named numeric vector: `ratea`, `rateb` (aggregate rates), `change`,
`rateeff`, `mixeff` (with `change = rateeff + mixeff`), `up` and `down`
(types whose rate rose or fell) and `simpson`: 2 when every type moved
one way and the aggregate the other, 1 when the aggregate change and the
rate effect differ in sign, 0 otherwise.

## Details

Each part's mix term is measured from the mean of the two aggregate
rates, so a part's term is negative when work moved toward a part that
clears below average; the centre changes how the mix effect is shared
among the parts, never its total.

## References

Kitagawa, E. M. (1955). Components of a difference between two rates.
*Journal of the American Statistical Association*, 50(272), 1168-1194.
[doi:10.1080/01621459.1955.10501299](https://doi.org/10.1080/01621459.1955.10501299)

## Examples

``` r
# both types clear better in B, yet the aggregate falls
clrsmp(clra = c(minor = 90, serious = 20), reca = c(100, 50),
       clrb = c(19, 60), recb = c(20, 140))
#> $parts
#>      part ratea     rateb    sharea shareb    rateeff     mixeff
#> 1   minor   0.9 0.9500000 0.6666667  0.125 0.01979167 -0.1687066
#> 2 serious   0.4 0.4285714 0.3333333  0.875 0.01726190 -0.1079303
#> 
#> $total
#>       ratea       rateb      change     rateeff      mixeff          up 
#>  0.73333333  0.49375000 -0.23958333  0.03705357 -0.27663690  2.00000000 
#>        down     simpson 
#>  0.00000000  2.00000000 
#> 
```
