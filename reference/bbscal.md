# Basic Space (blackbox) scaling

Recovers respondents' positions from a respondent-by-issue matrix of
responses by the decomposition of Poole (1998), \\X_0 = \Psi W' + J_n
c' + E_0\\ over the observed cells. Each dimension starts from the
sign-corrected respondent means (column polarities from the pairwise
correlations) and is refined by four alternating least-squares sweeps on
the residuals of the earlier dimensions; up to five joint sweeps follow,
stopping once the two half-step error sums of squares differ by less
than 0.01; the fitted \\\Psi W'\\ is then re-expressed by its singular
value decomposition, \\\Psi = U D^{1/2}\\, \\W = V D^{1/2}\\.
Respondents with fewer than `minscale` responses are not scaled. The
sign of each dimension is arbitrary.

## Usage

``` r
bbscal(X, ndim = 2, minscale = 8)
```

## Arguments

- X:

  A respondent-by-issue matrix of responses, `NA` for missing.

- ndim:

  The number of dimensions (at most one fewer than the issues and than
  the scaled respondents).

- minscale:

  The fewest responses a respondent needs to be scaled (capped at the
  number of issues).

## Value

A list of class `bbscal`: `ideal` (respondents by dimensions, `NA` rows
for those not scaled), `weights` (issues by dimensions), `eigen` and
`singular` (of the fitted \\\Psi W'\\), `fit` (\\1 - SSE/SS\\, SS about
the grand mean of the observed cells), `cmeans` (the issue intercepts
\\c\\) and `ndim`.

## References

Poole, K. T. (1998). Recovering a basic space from a set of issue
scales. *American Journal of Political Science*, 42(3), 954-993.
[doi:10.2307/2991737](https://doi.org/10.2307/2991737)

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md),
[`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md),
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
[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md),
[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)

## Examples

``` r
# 30 simulated respondents on eight issues, ten answers missing
set.seed(1)
X <- matrix(rnorm(60), 30, 2) %*% matrix(rnorm(16), 2, 8) +
  matrix(rnorm(240, 0, 0.3), 30, 8)
X[sample(240, 10)] <- NA
bbscal(X, ndim = 2)$fit
#> [1] 0.9787813
```
