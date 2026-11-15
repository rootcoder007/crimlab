# Wordfish: positions from word counts

The Poisson scaling model of Slapin and Proksch (2008), \\y\_{ij} \sim
\mathrm{Poisson}(\exp(\psi_i + \alpha_j + \beta_j \omega_i))\\, for
documents (judgments, briefs, platforms) by their word counts:
alternating Newton steps on the document positions \\\omega\\
(standardised after each sweep, the word parameters adjusted to keep the
linear predictor) and on the word parameters, from normal starts (seed
42, as [`set.seed()`](https://rdrr.io/r/base/Random.html) takes it),
until the positions change by a relative `tol`. The sign of the scale is
arbitrary.

## Usage

``` r
wrdfsh(dtm, maxit = 100, tol = 1e-06)
```

## Arguments

- dtm:

  A document-by-term matrix of counts.

- maxit:

  The most sweeps.

- tol:

  The convergence tolerance.

## Value

A list of class `wrdfsh`: `pos` (\\\omega\\), `beta`, `alpha`, `psi`,
`loglik` and `iter`.

## References

Slapin, J. B. and Proksch, S.-O. (2008). A scaling model for estimating
time-series party positions from texts. *American Journal of Political
Science*, 52(3), 705-722.
[doi:10.1111/j.1540-5907.2008.00338.x](https://doi.org/10.1111/j.1540-5907.2008.00338.x)

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
[`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md),
[`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md),
[`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md),
[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md)

## Examples

``` r
set.seed(1)
wrdfsh(matrix(rpois(600, 5), 20), maxit = 20)$pos
#>  [1]  1.15570659  0.25111445  0.01650556 -1.49106129  0.13088144 -1.07365520
#>  [7]  1.05617759 -0.83530154 -0.99936159  0.25353260 -0.31732379 -1.84826884
#> [13] -0.20769336  1.06437931  1.98755859  0.68871590  0.01462927 -0.99868738
#> [19]  1.10180015  0.05035155
```
