# Morisita's index of dispersion

Morisita's (1959) index of the counts in `n` equal quadrats,
`I = n (sum x^2 - S) / (S^2 - S)` with `S` the total: 1 for a random
pattern, above for clumping (`n` when all fall in one quadrat), below
for an even spread. The uniform and clumped indices are the values of
`I` at the `crit / 2` tails of chi-square on `n - 1` degrees of freedom,
and the standardised index of Smith-Gill (1975) maps `I` onto `[-1, 1]`
with them at `-0.5` and `0.5`; `chi = I (S - 1) + n - S` has the
chi-square distribution under randomness. These are the values of
`vegan::dispindmorisita()`.

## Usage

``` r
morsta(counts, crit = 0.05)
```

## Arguments

- counts:

  The counts: numbers from 0, at least two, adding up to more than 1 (a
  [`qdrtst()`](https://rootcoder007.github.io/crimlab/reference/qdrtst.md)
  count matrix will do).

- crit:

  The level of the critical indices, above 0 and below 1.

## Value

A named vector: `imor`, `mclu`, `muni`, `imst`, `chi` and `p` (its upper
chi-square tail).

## References

Morisita, M. (1959). Measuring of the dispersion of individuals and
analysis of the distributional patterns. *Memoirs of the Faculty of
Science, Kyushu University, Series E (Biology)*, 2, 215-235.

Smith-Gill, S. J. (1975). Cytophysiological basis of disruptive
pigmentary patterns in the leopard frog *Rana pipiens*. II. Wild type
and mutant cell-specific patterns. *Journal of Morphology*, 146(1),
35-54.
[doi:10.1002/jmor.1051460103](https://doi.org/10.1002/jmor.1051460103)

## Examples

``` r
morsta(c(3, 0, 7, 12, 1, 4, 0, 9))
#>         imor         mclu         muni         imst          chi            p 
#> 1.676190e+00 1.257508e+00 8.482820e-01 5.310481e-01 3.066667e+01 7.161935e-05 
```
