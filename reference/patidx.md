# Paternity index

The weight of the DNA evidence that an alleged father is the father of a
child, against an unrelated man, locus by locus and overall. With the
mother's genotype the index uses the allele the child must have received
from the father; without it the index is the parent-child
[`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md).
There is no allowance for mutation, so one inconsistent locus gives an
index of 0. See `proof("patidx")`.

## Usage

``` r
patidx(child, father, freq, mother = NULL)
```

## Arguments

- child, father:

  Genotypes: two-column matrices of allele labels, one row per locus, or
  two alleles for one locus.

- freq:

  The population frequencies of the child's alleles, same shape.

- mother:

  `NULL`, or the mother's genotype in the same form.

## Value

A list: `locus`, a data frame with `locus` and `idx` (`NaN` where the
mother's genotype rules her out as mother); `log10idx`; and `post`, the
probability of paternity at even prior odds.

## Examples

``` r
patidx(child = c(12, 14), father = c(14, 15), freq = c(0.1, 0.2),
       mother = c(12, 12))
#> $locus
#>   locus idx
#> 1     1 2.5
#> 
#> $log10idx
#> [1] 0.39794
#> 
#> $post
#> [1] 0.7142857
#> 
```
