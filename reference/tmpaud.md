# Disparity metrics over time and place

The disparate impact ratio, the parity gap, the Gini coefficient of the
groups' shares flagged and the bias amplification score of
[`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md),
for every city and period, and how far each city's ratio moves over the
periods (after Barman and Barman 2026).

## Usage

``` r
tmpaud(pred, grp = NULL, period = NULL, city = NULL, base = NULL)
```

## Arguments

- pred:

  The tool's decisions: 1 (or `TRUE`) for a case it flags, 0 (or
  `FALSE`) otherwise, one per case. A `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too: its elements fill the arguments by name.

- grp:

  Group labels: numbers, text or a factor, one per case, no missing
  values. Groups are ordered as the factor's levels, or sorted.

- period:

  The period of each case: labels (text, numbers or a factor), sorted as
  [`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md)
  sorts groups; ISO dates as text sort in time order.

- city:

  `NULL` (one place, labelled `"all"`), or the city of each case:
  labels.

- base:

  The base group: `NULL` for the group with the largest share flagged
  (the first in order on ties), its label, or its position from 1.

## Value

A list: `cells`, a data frame with one row per city and period that
occurs: `city`, `period`, `n`, `ok` (`TRUE` when audited), `di`, `dp`,
`gini` and `bas` (`NaN` when skipped); `cities`, one row per city:
`city`, `ncell` (cells audited), `di`, `dp`, `gini`, `bas` (their
means), `dimin`, `dimax`, `dirng` and `ngt1`; `sum`, `rng` (the largest
range) and `spread` (largest minus smallest mean ratio of the cities;
`NaN` with fewer than two); `base`.

## Details

A cell (one city in one period) is audited when it holds at least two
groups and the base group; otherwise it is reported and skipped. The
base group is the same in every cell: by default the group with the
largest share flagged over all the data. For each city, the means of the
four values over its audited cells (finite values only), the smallest
and largest ratio, their range, and how many periods have a ratio above

1.  A ratio that moves a lot from one period to the next says that one
    audit at one time describes that time only; a spread between cities
    says the same tool can treat groups differently in different places.
    Neither says why. For intervals on one cell, call
    [`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md)
    on it.

## References

Barman, P. K. and Barman, P. K. (2026). Spatial disparities in simulated
patrol detection: a reproducible generative-model audit across Baltimore
and Chicago.
[doi:10.48550/arXiv.2603.18987](https://doi.org/10.48550/arXiv.2603.18987)

## Examples

``` r
set.seed(2)
n <- 600
g <- sample(c("a", "b"), n, TRUE)
per <- sample(c("2024-1", "2024-2", "2024-3"), n, TRUE)
cty <- sample(c("north", "south"), n, TRUE)
p <- rbinom(n, 1, ifelse(g == "a", 0.4, ifelse(cty == "north", 0.3, 0.2)))
tmpaud(p, g, per, cty)$cities
#>    city ncell        di          dp       gini          bas     dimin     dimax
#> 1 north     3 0.7754133 -0.09992407 0.06547922 -0.008588778 0.6506747 0.9232143
#> 2 south     3 0.3603946 -0.28545866 0.23780500 -0.069788958 0.2545455 0.4545455
#>       dirng ngt1
#> 1 0.2725396    0
#> 2 0.2000000    0
```
