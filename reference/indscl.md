# INDSCAL: scaling with individual differences

The weighted Euclidean model of Carroll and Chang (1970): one stimulus
space shared by every individual, each individual weighting its
dimensions. The configuration starts from classical scaling of the mean
dissimilarities; then the weights of each individual and dimension (at
least 0.01) and the stimulus coordinates are updated in turn until the
configuration changes by a relative `tol`. The stress is the sum of
squared differences between each individual's dissimilarities and the
distances in its weighted space.

## Usage

``` r
indscl(D, ndim = 2, maxit = 300, tol = 1e-06)
```

## Arguments

- D:

  A list of square dissimilarity matrices of one size, one per
  individual (judge, panel, survey wave).

- ndim:

  The number of dimensions; `ndim` times the points scaled must be at
  most 67108864 (2^26).

- maxit:

  The most sweeps.

- tol:

  The convergence tolerance on the relative change of the configuration.

## Value

A list of class `indscl`: `config` (stimuli by dimensions), `weights`
(individuals by dimensions), `stress`, `iter`, `nind` and `nstim`.

## References

Carroll, J. D. and Chang, J.-J. (1970). Analysis of individual
differences in multidimensional scaling via an N-way generalization of
"Eckart-Young" decomposition. *Psychometrika*, 35(3), 283-319.
[doi:10.1007/BF02310791](https://doi.org/10.1007/BF02310791)

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
D1 <- as.matrix(dist(matrix(rnorm(20), 5)))
D2 <- as.matrix(dist(matrix(rnorm(20), 5)))
indscl(list(D1, D2), maxit = 30)$weights
#>      [,1]    [,2]
#> [1,]   NA 1979597
#> [2,]   NA 1929755
```
