# Least-squares metric unfolding

Unfolding of a respondent-by-stimulus matrix in the alternating
least-squares form of MLSMU6 (Poole 1984): respondents and stimuli take
gradient steps in turn against the row-centred dissimilarities, each
step scaled by the largest eigenvalue of the gradient's cross-product
and both clouds centred, until the stress
([`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md))
changes by a relative `tol`; of `nstart` normal starts (seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it) the lowest
stress is kept.

## Usage

``` r
mlsmu(D, ndim = 2, maxit = 200, tol = 1e-06, nstart = 5)
```

## Arguments

- D:

  A respondent-by-stimulus matrix of dissimilarities.

- ndim:

  The number of dimensions.

- maxit:

  The most sweeps per start.

- tol:

  The convergence tolerance on the relative change in stress.

- nstart:

  The random starts.

## Value

A list of class `mlsmu`: `xr`, `xs`, `stress`, `iter` and `conv`.

## References

Poole, K. T. (1984). Least squares metric, unidimensional unfolding.
*Psychometrika*, 49(3), 311-323.

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
D <- matrix(runif(120), 20, 6)
mlsmu(D, nstart = 1, maxit = 50)$stress
#> [1] 111.1204
```
