# Reliability of specimen matrices in a submerged or decomposed body

Ranks the specimen matrices of a case (vitreous humour, liver, femoral
blood) by how far their results can be trusted. Each matrix starts from
a reliability in a fresh body, `base`, which the user sets from the
literature for the analyte at hand; decomposition lowers every matrix by
`stage / smax`, at most by `scap`, and water dilutes the matrices marked
`dil` (blood, urine) by a further `days / dmax`, at most `dcap`. A
transparent ordering, not an estimate.

## Usage

``` r
mtxrel(
  base,
  dil,
  days = 0,
  stage = 0,
  dmax = 30,
  dcap = 0.8,
  smax = 10,
  scap = 0.6
)
```

## Arguments

- base:

  The reliability of each matrix in a fresh body: numbers from 0, named
  by matrix.

- dil:

  Whether water dilutes each matrix: `TRUE` or `FALSE`, one per matrix.

- days:

  The days the body was submerged: one number from 0.

- stage:

  The decomposition stage (0 fresh): one number from 0.

- dmax, smax:

  The days and the stage at which the loss reaches its rate's full
  scale: positive numbers.

- dcap, scap:

  The largest losses to dilution and to decomposition: numbers from 0 to
  1.

## Value

A data frame with one row per matrix, the most reliable first: `matrix`
(the names of `base`, or their numbers), `base`, `rel` the reliability
after the losses, `dil`, `rank`.

## References

Dinis-Oliveira, R. J., Carvalho, F., Duarte, J. A., Remiao, F., Marques,
A., Santos, A. and Magalhaes, T. (2010). Collection of biological
samples in forensic toxicology. *Toxicology Mechanisms and Methods*,
20(7), 363-414.
[doi:10.3109/15376516.2010.497976](https://doi.org/10.3109/15376516.2010.497976)

## Examples

``` r
base <- c(vitreous = 0.95, liver = 0.8, femoral = 0.55, heart = 0.35)
mtxrel(base, dil = c(FALSE, FALSE, TRUE, TRUE), days = 12, stage = 3)
#>     matrix base   rel   dil rank
#> 1 vitreous 0.95 0.665 FALSE    1
#> 2    liver 0.80 0.560 FALSE    2
#> 3  femoral 0.55 0.231  TRUE    3
#> 4    heart 0.35 0.147  TRUE    4
```
