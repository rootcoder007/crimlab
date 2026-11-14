# Pedigree likelihood ratio

The weight of DNA profiles for one family tree against another: is this
child the son of that man, are these remains the missing brother of the
people who reported him missing, are two people full or half siblings.
Each tree is a table of people and their parents; people who are not
typed (a dead parent, a missing relative) are summed over. The
likelihood of the profiles is computed under each tree, marker by
marker, and the ratio reported with its base-10 logarithm, the product
over independent markers.

## Usage

``` r
pedlr(ped, geno, freq, alt = NULL, mut = "none", rate = 0.001, range = 0.5)
```

## Arguments

- ped:

  The first tree: a data frame with columns `id`, `fid` and `mid` (the
  person, the father, the mother), one row per person; a founder has
  both parents missing, `""` or `0`. A parent must be a person of the
  table, no one may be both a father and a mother, and no one their own
  ancestor.

- geno:

  The profiles: a data frame with columns `id`, `locus`, `a1` and `a2`,
  one row per typed person and marker (labels; both alleles missing for
  a marker not typed).

- freq:

  The allele frequencies: a data frame with columns `locus`, `allele`
  and `freq` (above 0, at most 1, summing to at most 1 at each marker).
  Every allele of `geno` must be listed at its marker.

- alt:

  The second tree, laid out as `ped`; `NULL` (the default) for every
  typed person unrelated to every other.

- mut:

  The mutation model: `"none"`, `"equal"` or `"step"`.

- rate:

  The mutation rate per transmission: one number for both parents, or
  two, the maternal and the paternal rate.

- range:

  The ratio of the chance of a mutation of `d + 1` steps to that of `d`
  steps under `"step"`: one number strictly between 0 and 1.

## Value

An object of class `pedlr`, a list: `locus`, a data frame with `locus`,
`l1` and `l2` (the base-10 log likelihoods under `ped` and `alt`) and
`lr`, their ratio, for each marker of `freq` that someone is typed at;
`log10lr`, the base-10 logarithm of the product of the ratios; `lr`,
that product; and `mut`.

## Details

A founder (a person whose parents the table does not give) carries two
alleles drawn independently from the population frequencies
(Hardy-Weinberg, no coancestry); a child receives one of the two alleles
of each parent with probability 1/2, and that allele may change as it
passes, by the mutation model: `"none"`; `"equal"`, where an allele
mutates with probability `rate` to each other allele with equal chance;
or `"step"`, where it mutates with probability `rate` to the other
alleles in proportion to `range^d`, `d` the difference of the repeat
numbers (allele labels must then be whole numbers), so one-step
mutations are likeliest. These are the models of the Familias tradition
(Egeland et al. 2000; Simonsson and Mostad 2016); a mutation matrix is
used as it is, not made stationary, so a long chain of untyped people
can shift the allele frequencies slightly. Markers are taken as
independent (no linkage); the alleles not listed in `freq` at a marker
are taken together as one more allele (under `"equal"`, one more allele
to mutate to), which leaves the ratio unchanged without mutation; under
`"step"` the frequencies must sum to 1. Frequencies summing to within
1e-6 of 1 are divided by their sum.

The sum over the untyped people is exact: the alleles each person
received from the father and from the mother are eliminated one at a
time (Lauritzen and Sheehan 2003), which handles loops, inbreeding and
any number of generations; its cost grows with the number of untyped
people linked together and with the number of alleles, and a marker
whose elimination would need a table of more than 2^26 numbers is
refused. With `mut = "none"` and simple relationships the ratio equals
that of
[`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md)
and
[`patidx()`](https://rootcoder007.github.io/crimlab/reference/patidx.md);
mutation turns an inconsistent marker from a ratio of 0 into a small
positive one.

## References

Elston, R. C. and Stewart, J. (1971). A general model for the genetic
analysis of pedigree data. *Human Heredity*, 21(6), 523-542.
[doi:10.1159/000152448](https://doi.org/10.1159/000152448)

Egeland, T., Mostad, P. F., Mevag, B. and Stenersen, M. (2000). Beyond
traditional paternity and identification cases. *Forensic Science
International*, 110(1), 47-59.
[doi:10.1016/S0379-0738(00)00147-X](https://doi.org/10.1016/S0379-0738%2800%2900147-X)

Lauritzen, S. L. and Sheehan, N. A. (2003). Graphical models for genetic
analyses. *Statistical Science*, 18(4), 489-514.
[doi:10.1214/ss/1081443232](https://doi.org/10.1214/ss/1081443232)

Simonsson, I. and Mostad, P. (2016). Stationary mutation models.
*Forensic Science International: Genetics*, 23, 217-225.
[doi:10.1016/j.fsigen.2016.04.005](https://doi.org/10.1016/j.fsigen.2016.04.005)

## Examples

``` r
# is AF the father of C? The mother M is typed too.
trio <- data.frame(id = c("AF", "M", "C"), fid = c(NA, NA, "AF"),
                   mid = c(NA, NA, "M"))
alt <- data.frame(id = c("AF", "M", "C", "X"), fid = c(NA, NA, "X", NA),
                  mid = c(NA, NA, "M", NA))
geno <- data.frame(id = rep(c("AF", "M", "C"), 2),
                   locus = rep(c("D3", "vWA"), each = 3),
                   a1 = c(14, 15, 15, 16, 17, 17),
                   a2 = c(16, 17, 16, 18, 17, 19))
freq <- data.frame(locus = rep(c("D3", "vWA"), c(4, 4)),
                   allele = c(14, 15, 16, 17, 16, 17, 18, 19),
                   freq = c(0.1, 0.25, 0.25, 0.2, 0.2, 0.3, 0.2, 0.1))
pedlr(trio, geno, freq, alt)
#> pedlr: 2 markers, mutation model none
#>  locus       l1        l2 lr
#>     D3 -2.90309 -3.204120  2
#>    vWA     -Inf -3.142668  0
#> log10 LR -Inf; LR 0
# the child's vWA allele 19 cannot come from AF without a mutation
pedlr(trio, geno, freq, alt, mut = "equal", rate = 0.002)
#> pedlr: 2 markers, mutation model equal
#>  locus        l1        l2          lr
#>     D3 -2.904393 -3.204554 1.996002003
#>    vWA -5.444349 -3.141802 0.004982566
#> log10 LR -2.002386; LR 0.009945212
```
