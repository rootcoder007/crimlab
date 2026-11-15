# Ordinal IRT (ordinal factor model)

Each ordinal response is a thresholded latent normal, \\y^\*\_{ij} =
\lambda\_{j0} + \lambda_j' \phi_i + e\_{ij}\\, \\y\_{ij} = c\\ when
\\\gamma\_{j,c-1} \< y^\*\_{ij} \le \gamma\_{j,c}\\, with \\\gamma\_{j1}
= 0\\, \\\phi_i \sim N(0, I)\\, loadings \\\lambda_j \sim N(0, I /
L_0)\\ and a flat prior on the free cutpoints (Quinn 2004). Each sweep
updates the free cutpoints of every item by the Metropolis-Hastings step
of Cowles (1996) with \\y^\*\\ integrated out (step sizes tuned during
the burn-in only), then draws \\y^\*\\, the loadings and the ideal
points from their full conditionals. The chain starts from the leading
principal components of the mean-filled responses. Each draw is
reflected so the first item loads positively (rotated onto the first
draw with several dimensions). Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
ordirt(Y, ndim = 1, nsamp = 500, burn = 100, seed = 42, lprec = 0)
```

## Arguments

- Y:

  A respondent-by-item matrix of ordinal responses (any ordered numeric
  coding), `NA` for missing.

- ndim:

  The number of dimensions.

- nsamp:

  The draws kept.

- burn:

  The sweeps discarded first.

- seed:

  The seed.

- lprec:

  The prior precision \\L_0\\ of the loadings and intercepts (0, flat).

## Value

A list of class `ordirt`: `ideal` and `sd`, `disc` (items by
dimensions), `icept`, `cuts` (one vector per item, the first 0),
`accept` (cutpoint acceptance after the burn-in) and `nsamp`.

## References

Quinn, K. M. (2004). Bayesian factor analysis for mixed ordinal and
continuous responses. *Political Analysis*, 12(4), 338-353.

Cowles, M. K. (1996). Accelerating Monte Carlo Markov chain convergence
for cumulative-link generalized linear models. *Statistics and
Computing*, 6(2), 101-111.

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md),
[`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md),
[`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md),
[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md),
[`bunfld()`](https://rootcoder007.github.io/crimlab/reference/bunfld.md),
[`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md),
[`cmdsf()`](https://rootcoder007.github.io/crimlab/reference/cmdsf.md),
[`cutlns()`](https://rootcoder007.github.io/crimlab/reference/cutlns.md),
[`dblctr()`](https://rootcoder007.github.io/crimlab/reference/dblctr.md),
[`dwnom()`](https://rootcoder007.github.io/crimlab/reference/dwnom.md),
[`dynirt()`](https://rootcoder007.github.io/crimlab/reference/dynirt.md),
[`indscl()`](https://rootcoder007.github.io/crimlab/reference/indscl.md),
[`ipsrec()`](https://rootcoder007.github.io/crimlab/reference/ipsrec.md),
[`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md),
[`irtlik()`](https://rootcoder007.github.io/crimlab/reference/irtlik.md),
[`irtpst()`](https://rootcoder007.github.io/crimlab/reference/irtpst.md),
[`mdsfit()`](https://rootcoder007.github.io/crimlab/reference/mdsfit.md),
[`mlsmu()`](https://rootcoder007.github.io/crimlab/reference/mlsmu.md),
[`nmdsf()`](https://rootcoder007.github.io/crimlab/reference/nmdsf.md),
[`nombt()`](https://rootcoder007.github.io/crimlab/reference/nombt.md),
[`nomll()`](https://rootcoder007.github.io/crimlab/reference/nomll.md),
[`nompr()`](https://rootcoder007.github.io/crimlab/reference/nompr.md),
[`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md),
[`nrmvec()`](https://rootcoder007.github.io/crimlab/reference/nrmvec.md),
[`optcls()`](https://rootcoder007.github.io/crimlab/reference/optcls.md),
[`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md),
[`procrs()`](https://rootcoder007.github.io/crimlab/reference/procrs.md),
[`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md),
[`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md),
[`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md),
[`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md),
[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md),
[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)

## Examples

``` r
set.seed(1)
th <- rnorm(60)
Y <- matrix(cut(outer(th, runif(6, 0.8, 1.5)) + rnorm(360),
                c(-Inf, -0.5, 0.5, Inf), labels = FALSE), 60, 6)
f <- ordirt(Y, nsamp = 50, burn = 50)
cor(f$ideal[, 1], th)
#> [1] 0.8995402
```
