# Genotype match probability

The chance that an unrelated person from the same population shares a
DNA profile, locus by locus and overall. With `theta = 0` this is the
product rule (`p^2` for a homozygote, `2pq` for a heterozygote). A
positive `theta` applies the Balding-Nichols correction for shared
ancestry. That never makes a homozygote match rarer than the product
rule says; it can make a heterozygote match rarer when both alleles are
common. See `proof("gtmp")`.

## Usage

``` r
gtmp(p, q = NULL, hom = FALSE, theta = 0)
```

## Arguments

- p:

  Frequency of the first allele at each locus: numbers above 0 and at
  most 1, optionally named by locus, as a numeric vector, a 1-d table, a
  list of single numbers or a one-row data frame.

- q:

  Frequency of the second allele at each locus; not read where `hom` is
  `TRUE`. Named inputs are matched to `p` by locus name, unnamed ones by
  position (only names protect against a mismatch).

- hom:

  `TRUE` where the genotype is homozygous: one value, or one per locus,
  matched to `p` by name when both are named.

- theta:

  The coancestry coefficient: one number from 0 up to, not including, 1.

## Value

A list: `locus`, a data frame with `locus`, `hom` and `mp`; `log10mp`,
the base-10 logarithm of the profile match probability; and `log10lr`,
the same with the sign changed, the weight of a single-source match
against an unrelated alternative.

## References

Balding, D. J. and Nichols, R. A. (1994). DNA profile match probability
calculation. *Forensic Science International*, 64(2-3), 125-140.
[doi:10.1016/0379-0738(94)90222-4](https://doi.org/10.1016/0379-0738%2894%2990222-4)

## Examples

``` r
gtmp(p = c(D3 = 0.1, vWA = 0.2), q = c(0.3, NA), hom = c(FALSE, TRUE),
     theta = 0.03)
#> $locus
#>   locus   hom         mp
#> 1    D3 FALSE 0.07467851
#> 2   vWA  TRUE 0.06607071
#> 
#> $log10mp
#> [1] -2.306795
#> 
#> $log10lr
#> [1] 2.306795
#> 
```
