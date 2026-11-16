# SMACOF stress minimisation

Metric scaling by majorisation (de Leeuw 1977): each step is the Guttman
transform \\X \leftarrow V^{+} B(X) X\\ with \\V = \mathrm{diag}(W 1) -
W\\, until the raw stress \\\frac{1}{2} \sum w\_{ij} (d\_{ij} -
\delta\_{ij})^2\\ changes by less than `tol`.

## Usage

``` r
smacf(D, ndim = 2, maxit = 300, tol = 1e-06, w = NULL, init = NULL)
```

## Arguments

- D:

  A symmetric matrix of dissimilarities.

- ndim:

  The number of dimensions; `ndim` times the points scaled must be at
  most 67108864 (2^26).

- maxit:

  The most steps.

- tol:

  The convergence tolerance on the change in stress.

- w:

  `NULL` for equal weights, or a matrix of weights.

- init:

  `NULL` for a random start (normal, from seed 42 as
  [`set.seed()`](https://rdrr.io/r/base/Random.html) takes it), or a
  starting configuration.

## Value

A list of class `smacf`: `coords`, `stress`, `iter` and `conv`.

## References

de Leeuw, J. (1977). Applications of convex analysis to multidimensional
scaling. In J. R. Barra et al. (eds.), *Recent Developments in
Statistics*, 133-145. North-Holland.

Borg, I. and Groenen, P. J. F. (2005). *Modern Multidimensional
Scaling*, 2nd ed. Springer.
[doi:10.1007/0-387-28981-X](https://doi.org/10.1007/0-387-28981-X)

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
[`ordirt()`](https://rootcoder007.github.io/crimlab/reference/ordirt.md),
[`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md),
[`procrs()`](https://rootcoder007.github.io/crimlab/reference/procrs.md),
[`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md),
[`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md),
[`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md),
[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md),
[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)

## Examples

``` r
set.seed(1)
D <- as.matrix(dist(matrix(rnorm(40), 10)))
smacf(D)$stress
#> [1] 6.626513
```
