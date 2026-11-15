# Trimmed inverse-probability weights (C backend)

The inverse-probability-of-treatment weights \\1/e\\ for the treated and
\\1/(1-e)\\ for the untreated, with the propensity score clamped into
`[trimlo, trimhi]` first: an untrimmed score near 0 or 1 gives one
observation a weight large enough to dominate the estimate. The weights
reweight a treated and an untreated group so that each resembles the
whole sample on the covariates behind the propensity score (Rosenbaum
and Rubin, 1983).

## Usage

``` r
cipw(treat, ps, trimlo = 0.01, trimhi = 0.99)
```

## Arguments

- treat:

  Treatment indicator, 0/1 or `TRUE`/`FALSE`; `NA` allowed.

- ps:

  Propensity scores in \\\[0, 1\]\\, the length of `treat`; `NA`
  allowed.

- trimlo, trimhi:

  The clamp bounds, `0 < trimlo < trimhi < 1`.

## Value

The weights, the length of `treat`.

## Details

With the clamped score \\\tilde e_i = \min(\max(e_i, l), u)\\, the
weight is \\1 / \tilde e_i\\ when \\T_i = 1\\ and \\1 / (1 - \tilde
e_i)\\ when \\T_i = 0\\; with the defaults none exceeds 100. A treatment
coded other than 0/1, a propensity outside \\\[0, 1\]\\ and bounds
without \\0 \< l \< u \< 1\\ are refused. A missing treatment or
propensity gives a missing weight for that unit only. The weights are
not normalised.

## References

Rosenbaum PR, Rubin DB (1983). The central role of the propensity score
in observational studies for causal effects. *Biometrika* 70(1), 41–55.
[doi:10.1093/biomet/70.1.41](https://doi.org/10.1093/biomet/70.1.41)

## See also

[`cwtd()`](https://rootcoder007.github.io/crimlab/reference/cwtd.md).

## Examples

``` r
cipw(c(1, 0, 1, 0), c(0.5, 0.25, 0.02, 0.9))
#> [1]  2.000000  1.333333 50.000000 10.000000
cipw(c(1, 0, 1, 0), c(0.5, 0.25, 0.02, 0.9), trimlo = 0.1)
#> [1]  2.000000  1.333333 10.000000 10.000000
cipw(c(TRUE, FALSE), c(0.4, 0.4))
#> [1] 2.500000 1.666667
try(cipw(c(1, 2), c(0.5, 0.5)))
#> Error : cipw: `treat` must be coded 0/1.
```
