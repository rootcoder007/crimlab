# Risk difference when group membership is misrecorded

When the group a person is recorded in (an officer's perception of race,
say) differs at random from the true one, with the same error rates
whatever the outcome, the measured groups are mixtures of the true ones.
The measured risk difference is then the true one times the factor
`(se + sp - 1) p (1 - p) / (d1 d0)`, where `d1` and `d0` are the shares
recorded in each group. For any sensitivity and specificity it is never
larger in size than the true difference, so a disparity in recorded data
understates the true one; with `se + sp >= 1` it never has the opposite
sign, while with `se + sp < 1` a difference other than 0 has it. See
`proof("misrd")`.

## Usage

``` r
misrd(r1, r0, prev, se, sp)
```

## Arguments

- r1:

  The true risk in group 1: one number from 0 to 1.

- r0:

  The true risk in group 0: one number from 0 to 1.

- prev:

  The true share in group 1, strictly between 0 and 1.

- se:

  Sensitivity: the chance a member of group 1 is recorded in group 1,
  from 0 to 1.

- sp:

  Specificity: the chance a member of group 0 is recorded in group 0,
  from 0 to 1.

## Value

A named numeric vector: `true` (`r1 - r0`), `measured` (the difference
between the recorded groups), `factor` (`measured / true` whenever
`true` is not 0), `d1` and `d0`.

## References

Bross, I. (1954). Misclassification in 2 x 2 tables. *Biometrics*,
10(4), 478-486. [doi:10.2307/3001619](https://doi.org/10.2307/3001619)

## Examples

``` r
misrd(r1 = 0.3, r0 = 0.1, prev = 0.4, se = 0.8, sp = 0.9)
#>      true  measured    factor        d1        d0 
#> 0.2000000 0.1426146 0.7130730 0.3800000 0.6200000 
```
