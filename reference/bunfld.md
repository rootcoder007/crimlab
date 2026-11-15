# Bayesian unfolding

The lognormal model of
[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md)
for a respondent-by-stimulus matrix: \\\log \delta\_{ij} \sim N(\log
\lVert x_i - z_j \rVert, 1 / \tau)\\, \\N(0, 10^2)\\ priors on every
coordinate and \\\tau \sim U(0, 10)\\ (Bakker and Poole 2013). Feeling
thermometers \\T\\ on 0-100 enter as \\\delta = (100 - T) / 50\\. Given
the stimuli the respondents are independent, and the other way round, so
each block is slice sampled at once, dimension by dimension; draws are
aligned on the stimuli. The chain starts from the singular value
decomposition of the centred matrix. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
bunfld(D, ndim = 2, nsamp = 1000, burn = 200, seed = 42)
```

## Arguments

- D:

  A respondent-by-stimulus matrix of dissimilarities.

- ndim:

  The number of dimensions.

- nsamp:

  The draws kept.

- burn:

  The sweeps discarded first.

- seed:

  The seed.

## Value

A list of class `bunfld`: `stim` (aligned posterior means), `sd`,
`ideal`, `dist` (posterior mean respondent-stimulus distances), `sigma`,
`tau` and `nsamp`.

## References

Bakker, R. and Poole, K. T. (2013). Bayesian metric multidimensional
scaling. *Political Analysis*, 21(1), 125-140.
[doi:10.1093/pan/mps039](https://doi.org/10.1093/pan/mps039)

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md),
[`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md),
[`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md),
[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md),
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
[`ordirt()`](https://rootcoder007.github.io/crimlab/reference/ordirt.md),
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
S <- matrix(runif(8, -0.5, 0.5), 4)
R <- matrix(runif(60, -0.5, 0.5), 30)
D <- sqrt(outer(R[, 1], S[, 1], "-")^2 + outer(R[, 2], S[, 2], "-")^2)
f <- bunfld(D, nsamp = 50, burn = 50)
cor(as.numeric(f$dist), as.numeric(D))
#> [1] 0.9766247
```
