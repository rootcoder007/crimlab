# Population stability index and Jensen-Shannon divergence (C backend)

Two effect sizes of how far the distribution of a numeric column has
moved from a reference sample. Both samples are binned on the
reference's quantiles and the two sets of bin shares compared. Unlike a
p-value they do not grow with the sample size, so they say how large a
shift is rather than whether it is detectable.

## Usage

``` r
drftpsi(x, y, bins = 10L, eps = 1e-06)
```

## Arguments

- x, y:

  Numeric vectors, the reference and the new sample; logical vectors are
  read as 0/1.

- bins:

  The number of quantile bins of the reference, one whole number from 2
  to 1,000,000; fewer are used when quantiles tie.

- eps:

  The floor of an empty bin in the index, in (0, 1). It is part of the
  statistic: compare runs at one `eps`.

## Value

A named numeric vector of class `drftpsi`: `psi` and `jsdiv`, with the
inputs in the attribute `lakner`.

## Details

Binning: values are rounded to 12 significant digits (so a formatting
difference between two releases does not move tied values across an
edge), the `bins + 1` type-7 quantiles of `x` are computed, duplicates
removed, and the interior quantiles become right-closed cut points with
open ends. When the reference has fewer than three distinct quantiles,
every distinct quantile is a cut point instead.

The population stability index is \\\sum_i (p_i - q_i)\log(p_i/q_i)\\
over the reference shares \\p_i\\ and the new shares \\q_i\\: the
symmetrised Kullback-Leibler divergence. The conventional reading, from
credit-risk monitoring where it originates, is below 0.1 stable, 0.1 to
0.25 worth a look, above 0.25 a material shift; these bands are
heuristics with no calibrated null distribution, and on a few hundred
rows binning noise alone can exceed 0.25. The Jensen-Shannon divergence
\\\tfrac12 KL(p\\m) + \tfrac12 KL(q\\m)\\, \\m = (p+q)/2\\, is bounded
by \\\log 2\\ nats, so it compares across columns with different numbers
of bins, and is finite when a bin is empty on one side. Empty bins are
floored at `eps` for the index only.

## References

Wu D, Olson DL (2010). Enterprise risk management: coping with model
risk in a large bank. *Journal of the Operational Research Society*
61(2), 179–190.
[doi:10.1057/jors.2008.144](https://doi.org/10.1057/jors.2008.144)

Lin J (1991). Divergence measures based on the Shannon entropy. *IEEE
Transactions on Information Theory* 37(1), 145–151.
[doi:10.1109/18.61115](https://doi.org/10.1109/18.61115)

## See also

[`drftks()`](https://rootcoder007.github.io/crimlab/reference/drftks.md).

## Examples

``` r
set.seed(2)
ref <- stats::rnorm(500)
drftpsi(ref, stats::rnorm(500))
#> drftpsi: psi = 0.022072, jsdiv = 0.002753 
#> inputs: nx = 500, ny = 500, bins = 10 
drftpsi(ref, stats::rnorm(500, mean = 1))
#> drftpsi: psi = 0.9543, jsdiv = 0.1053 
#> inputs: nx = 500, ny = 500, bins = 10 
drftpsi(ref, ref)
#> drftpsi: psi = 0, jsdiv = 0 
#> inputs: nx = 500, ny = 500, bins = 10 
drftpsi(c(1, 1, 1), c(9, 9, 9))[["jsdiv"]] <= log(2)
#> [1] TRUE
drftpsi(ref, stats::rnorm(500, mean = 1), bins = 4, eps = 0.01)
#> drftpsi: psi = 0.69315, jsdiv = 0.08235 
#> inputs: nx = 500, ny = 500, bins = 4 
try(drftpsi(ref, ref, bins = 1))
#> Error : drftpsi: `bins` must be one whole number of at least 2.
```
