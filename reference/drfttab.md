# Drift of a whole table between two releases

Compares a reference table and a current one column by column, with the
tests of this family: a numeric (or date) column by the two-sample
Kolmogorov-Smirnov test
([`drftks()`](https://rootcoder007.github.io/crimlab/reference/drftks.md))
and the population stability index
([`drftpsi()`](https://rootcoder007.github.io/crimlab/reference/drftpsi.md)),
a text, factor or logical column by the chi-square test of homogeneity
([`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md)).
`drftcal()` then says how often these screens fire on a table that has
not changed at all, which is what a flag has to be read against.

## Usage

``` r
drfttab(
  ref,
  cur,
  alpha = 0.01,
  psithr = 0.25,
  psiminn = 100L,
  bins = 10L,
  seed = 1L
)

drftcal(
  data,
  n = 50L,
  alpha = 0.01,
  seed = 1L,
  psithr = 0.25,
  psiminn = 100L,
  bins = 10L
)
```

## Arguments

- ref, cur:

  The reference and the current table: data frames. Only the columns
  named in both are compared.

- alpha:

  The level of each column's test, in (0, 1).

- psithr:

  The population stability index above which a numeric column is
  flagged; 0.25 is the conventional "material shift".

- psiminn:

  The fewest values each sample needs for `psithr` to count.

- bins:

  The number of quantile bins of the reference, one whole number from 2
  to 1,000,000; fewer are used when quantiles tie.

- seed:

  The seed of the Monte Carlo p-value
  [`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md)
  uses for small expected counts, and in `drftcal()` of the random
  splits too (there a whole number from 0:
  [`cunif()`](https://rootcoder007.github.io/crimlab/reference/cunif.md)
  draws a seed for each split from it, and a split takes the rows of the
  smaller half of the uniforms drawn at its seed, so a seed gives the
  same splits in every language of the package).

- data:

  One table, known not to have drifted (one release).

- n:

  The number of random half-splits, from 1 to 1,000,000.

## Value

`drfttab()` an object of class `drfttab`, a data frame with a row per
shared column: `column`, `type` (`"numeric"` or `"label"`), `test`,
`statistic`, `p`, `psi` (`NA` for a label column), `nref`, `ncur`,
`flag` and `note`; the attributes `only` (the columns in one table
alone), `alpha`, `psithr` and `psiminn`. `drftcal()` an object of class
`drftcal`, a list: `columns`, a data frame of `column` and `far` (the
share of splits on which it was flagged); `anyflag`, the share of splits
with at least one flag; `alpha`; `n`; and `alphafam`, `alpha` divided by
the number of columns, the level per column that holds the chance of any
false flag at `alpha` (Bonferroni).

## Details

`drfttab()` screens the columns the two tables share. A column is
flagged when its test rejects at `alpha`, or when its population
stability index is above `psithr` and both samples have at least
`psiminn` values (on fewer, binning noise alone crosses the usual
bands). A column that cannot be tested (all missing on one side, one
category, of a type that is neither numbers nor labels) has `flag` `NA`,
and `note` says why. Missing values are left out of each test; `nref`
and `ncur` are the counts that entered it.

Many columns mean many tests: at `alpha` 0.01, one column in a hundred
is flagged on unchanged data by the tests alone. `drftcal()` measures
this for the table at hand: it splits `data` at random into two halves
`n` times, runs `drfttab()` on each pair, and reports the share of
splits on which each column, and any column, was flagged. Each half has
half the rows, so the rates describe a comparison at that size.

Both functions exist in every language of the package. One branch is R's
alone: a column of labels with an expected count below 5 gets its
p-value from a Monte Carlo draw of R's generator
([`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md));
the other languages leave such a column untested (`flag` missing) and
say so in `note`.

Errors: `ref`, `cur` or `data` not a data frame; no shared column;
`alpha` outside (0, 1); `psithr` not a positive number; `psiminn`,
`bins` or `n` not a whole number of at least 1 (`bins` at least 2, `n`
at most 1,000,000); fewer than four rows in `data`; a `seed` that is not
one whole number (`drftcal()`: from 0).

## References

Dunn, O. J. (1961). Multiple comparisons among means. *Journal of the
American Statistical Association*, 56(293), 52-64.
[doi:10.1080/01621459.1961.10482090](https://doi.org/10.1080/01621459.1961.10482090)

## See also

[`drftks()`](https://rootcoder007.github.io/crimlab/reference/drftks.md),
[`drftpsi()`](https://rootcoder007.github.io/crimlab/reference/drftpsi.md),
[`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md)
for one column.

## Examples

``` r
set.seed(1)
ref <- data.frame(age = rnorm(300, 35, 9),
                  charge = sample(c("assault", "theft", "drugs"), 300, TRUE))
cur <- data.frame(age = rnorm(300, 39, 9),
                  charge = sample(c("assault", "theft", "drugs"), 300, TRUE))
drfttab(ref, cur)
#> drfttab: 2 columns, 1 flagged, 0 not testable (alpha 0.01, PSI above 0.25)
#>  column    type        test statistic        p   psi nref ncur  flag
#>     age numeric          ks     0.200 1.23e-05 0.309  300  300  TRUE
#>  charge   label homogeneity     2.029 3.63e-01    NA  300  300 FALSE
# how often the screens fire when nothing changed
drftcal(ref, n = 20)
#> drftcal: 20 half-splits of unchanged data at alpha 0.01
#> any column flagged on 5% of splits; per-column level for a family-wise 0.01: 0.005
#>  column  far
#>  charge 0.05
try(drfttab(ref, data.frame(other = 1)))
#> Error : drfttab: `ref` and `cur` share no column.
```
