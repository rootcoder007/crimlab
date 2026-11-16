# Randomisation lists for a trial

Allocates participants to arms by simple, permuted-block or stratified
block randomisation from a seed, so that the list can be made again
exactly from the seed and the settings, in R or on the command line.

## Usage

``` r
rndalc(
  n = NULL,
  arms = c("A", "B"),
  ratio = rep(1, length(arms)),
  block = c(2, 3) * sum(ratio),
  strata = NULL,
  seed = NULL,
  prefix = "P"
)

# S3 method for class 'rndalc'
print(x, ...)
```

## Arguments

- n:

  The number of participants: a whole number from 1 to 10000000. With
  `strata` it may be left out (it is their total).

- arms:

  The names of the arms: two or more, each once.

- ratio:

  The allocation ratio: one whole number from 1 per arm.

- block:

  The block sizes drawn from: whole numbers, each a multiple of
  `sum(ratio)`; 0 for simple randomisation.

- strata:

  `NULL`, or the strata and their places: a named vector of whole
  numbers from 0.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647.

- prefix:

  The start of the participant ids: one string.

- x:

  An object of class `rndalc`.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `rndalc`, a list: `alloc`, a data frame in enrolment
order with `id`, `stratum` (with `strata`), `block` (the block's number
within its stratum; `NA` for simple randomisation) and `arm`; `tab`, the
participants of each arm by stratum (a matrix, one row `all` without
strata); and `seed`, the seed used.

## Details

With `block = 0` each participant is allocated alone (simple
randomisation), to arm `a` with chance `ratio[a] / sum(ratio)`.
Otherwise blocks are made one after another, each of a size drawn with
equal chances from `block` (every size a multiple of `sum(ratio)`),
holding the arms in the ratio in an order drawn uniformly (a
Fisher-Yates shuffle). Block sizes drawn at random keep the next
allocation hard to guess at a site (Schulz and Grimes 2002). The list
stops at the places asked for, so the last block of a stratum may be
incomplete. With `strata`, each stratum has its own sequence of blocks
(Rosenberger and Lachin 2016, ch. 3).

The draws come from the package's own generator (xorshift64\* seeded
through splitmix64; Vigna 2016), the same on every platform; the seed is
returned with the list. Sealing the list (a commitment signed and
timestamped before enrolment) is left to other tools.

## References

Rosenberger, W. F. and Lachin, J. M. (2016). *Randomization in Clinical
Trials: Theory and Practice*, 2nd ed. Wiley.
[doi:10.1002/9781118742112](https://doi.org/10.1002/9781118742112)

Schulz, K. F. and Grimes, D. A. (2002). Generation of allocation
sequences in randomised trials: chance, not choice. *The Lancet*,
359(9305), 515-519.
[doi:10.1016/S0140-6736(02)07683-3](https://doi.org/10.1016/S0140-6736%2802%2907683-3)

Vigna, S. (2016). An experimental exploration of Marsaglia's xorshift
generators, scrambled. *ACM Transactions on Mathematical Software*,
42(4), 30. [doi:10.1145/2845077](https://doi.org/10.1145/2845077)

## Examples

``` r
a <- rndalc(40, seed = 12)
a
#> rndalc: 40 participants, seed 12
#>      A  B
#> all 20 20
head(a$alloc)
#>    id block arm
#> 1 P01     1   A
#> 2 P02     1   B
#> 3 P03     1   A
#> 4 P04     1   A
#> 5 P05     1   B
#> 6 P06     1   B
s <- rndalc(arms = c("drug", "placebo"), ratio = c(2, 1),
            strata = c(Toronto = 30, Ottawa = 24), seed = 7)
s$tab
#>         drug placebo
#> Toronto   21       9
#> Ottawa    16       8
rndalc(20, block = 0, seed = 1)$tab
#>      A B
#> all 13 7
```
