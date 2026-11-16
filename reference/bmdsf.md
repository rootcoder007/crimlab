# Bayesian metric multidimensional scaling

Each observed dissimilarity is lognormal about the distance between the
two objects' coordinates, \\\log \delta\_{ij} \sim N(\log d\_{ij}, 1 /
\tau)\\, with \\x\_{ik} \sim N(0, 10^2)\\ and \\\tau \sim U(0, 10)\\
(Bakker and Poole 2013). Coordinates are slice sampled one at a time
(Neal 2003), \\\tau\\ drawn from its truncated gamma; the chain starts
from classical scaling plus small normal noise. The likelihood is
unchanged by translation, rotation and reflection, so every draw is
aligned (translation and rotation) onto the posterior mean. Zero,
negative or missing dissimilarities are left out. Draws come from
`seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
bmdsf(D, ndim = 2, nsamp = 1000, burn = 200, sigma = 1, seed = 42)
```

## Arguments

- D:

  A square matrix of dissimilarities.

- ndim:

  The number of dimensions; `ndim` times the points scaled must be at
  most 67108864 (2^26).

- nsamp:

  The draws kept; `nsamp` times the values kept per draw must be at most
  67108864 (2^26).

- burn:

  The sweeps discarded first; `burn` plus `nsamp` must be at most
  2147483647.

- sigma:

  The starting \\\sigma = \tau^{-1/2}\\.

- seed:

  The seed.

## Value

A list of class `bmdsf`: `pos` (aligned posterior means), `sd`, `dist`
(posterior mean distances), `sigma`, `tau`, `draws` (draws by objects by
dimensions) and `nsamp`.

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
D <- as.matrix(dist(matrix(rnorm(16), 8)))
f <- bmdsf(D, nsamp = 50, burn = 50)
cor(f$dist[upper.tri(D)], D[upper.tri(D)])
#> [1] 0.9757443
```
