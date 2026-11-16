# Simulated records with a known disparity

Incidents in segregated areas, whose true occurrence is the same for
every group but whose recording, and the risk an algorithm predicts for
them, differ by group by a set amount: data on which an audit's answer
is known (after Barman and Barman 2026). The first group is the
reference.

## Usage

``` r
fairsim(
  n = 2000,
  groups = c("A", "B"),
  props = NULL,
  nareas = 20,
  base = 0.3,
  bias = 0.5,
  seed = NULL
)
```

## Arguments

- n:

  The number of records: a whole number from 1.

- groups:

  The group labels, at least two, the reference first.

- props:

  `NULL` for equal shares, or one weight from 0 per group (not all 0).

- nareas:

  The number of areas: a whole number, at least the number of groups.

- base:

  The true rate, and the reference group's recorded rate: one number
  from 0 to 1.

- bias:

  The disparity injected: one number from -1 to 1.

- seed:

  `NULL` to draw the seed from R's generator, or one whole number from
  0; the records come from crimlab's own generator, the same on every
  platform.

## Value

A data frame with one row per record: `area` (a number from 1), `group`,
`truth`, `det` and `risk`.

## Details

Each record draws its group by `props`, then an area among those its
group lives in (with `g` groups, areas 1, `g + 1`, `2 g + 1`, ... hold
the first group only, areas 2, `g + 2`, ... the second, and so on), then
`truth`, 1 with probability `base` whatever the group; `det`, the
recorded outcome, 1 with probability `base` in the reference group and
`base (1 - bias)` in the others (clamped to 0 and 1), so that the other
groups' incidents are recorded less when `bias` is positive; and `risk`,
normal with mean 250 in the reference group and `250 + 100 bias` in the
others and standard deviation 40, clamped to 0 and 500. So
[`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md)
on `det` should find a ratio near `1 - bias`,
[`scrdsp()`](https://rootcoder007.github.io/crimlab/reference/scrdsp.md)
on `risk` a gap near `100 bias`, and
[`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md)
the other groups' areas ranked as riskier than their recorded outcomes.

## References

Barman, P. K. and Barman, P. K. (2026). Spatial disparities in simulated
patrol detection: a reproducible generative-model audit across Baltimore
and Chicago.
[doi:10.48550/arXiv.2603.18987](https://doi.org/10.48550/arXiv.2603.18987)

## Examples

``` r
d <- fairsim(2000, seed = 1)
head(d)
#>   area group truth det     risk
#> 1   10     B     0   0 379.1753
#> 2    2     B     0   1 324.0614
#> 3   10     B     0   1 280.4492
#> 4    1     A     0   0 312.1302
#> 5    6     B     1   0 326.8925
#> 6    8     B     0   1 304.4651
# the injected ratio is 1 - 0.5
fairdf(d$det, d$group, base = "A")$tab[, c("group", "di", "dilo", "dihi")]
#>   group        di      dilo      dihi
#> 1     A 1.0000000       NaN       NaN
#> 2     B 0.5751675 0.4849938 0.6821069
```
