# Kinship index

The weight of two DNA profiles for a stated relationship against no
relationship, locus by locus and overall. Loci are taken as independent;
there is no allowance for mutation or coancestry. See `proof("kinidx")`.

## Usage

``` r
kinidx(g1, g2, freq, rel = "parent")
```

## Arguments

- g1, g2:

  The two genotypes: two-column matrices of allele labels, one row per
  locus (row names name the loci), or two alleles for one locus.

- freq:

  The population frequencies of the alleles of `g1`, same shape.

- rel:

  The relationship: `"parent"`, `"sibling"`, `"half"` (half siblings,
  grandparent or uncle) or `"cousin"`, or the three probabilities of
  sharing 0, 1 and 2 alleles by descent. Names are matched exactly, in
  lower case.

## Value

A list: `locus`, a data frame with `locus` and `idx`; `log10idx`, the
base-10 logarithm of the combined index; and `post`, the probability of
the relationship at even prior odds.

## Examples

``` r
kinidx(g1 = c(12, 14), g2 = c(12, 14), freq = c(0.1, 0.2), rel = "sibling")
#> $locus
#>   locus   idx
#> 1     1 8.375
#> 
#> $log10idx
#> [1] 0.9229848
#> 
#> $post
#> [1] 0.8933333
#> 
```
