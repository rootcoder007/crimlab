# Aldrich-McKelvey scaling

Recovers the positions of stimuli (parties, candidates, courts) on one
scale from respondents' placements of them, each respondent allowed an
own shift and stretch of the scale: \\z\_{ij} = a_i + b_i \hat z_j +
\epsilon\_{ij}\\. The positions are the least-squares solution of
Aldrich and McKelvey (1977): with \\X_i = \[1, z_i\]\\ and \\P_i\\ the
projection onto its columns, \\\hat z\\ is the eigenvector of \\\sum_i
(I - P_i)\\ with the smallest eigenvalue orthogonal to the constant.
Respondents with a missing placement, or who place every stimulus at the
same point, are left out of that sum. The positions are standardised to
mean 0 and standard deviation 1, the first stimulus on the left. The
respondent intercepts and slopes are each respondent's least-squares
regression of the placements on those positions (defined with two
placements or more; a slope of zero is set to 1e-10).

## Usage

``` r
amscal(Z)
```

## Arguments

- Z:

  A respondent-by-stimulus matrix of placements, `NA` for a missing one.

## Value

A list of class `amscal`: `zhat` (stimulus positions), `alpha` and
`beta` (each respondent's intercept and slope) and `weights`
(\\\|b_i\|\\ scaled to sum to the number of respondents).

## References

Aldrich, J. H. and McKelvey, R. D. (1977). A method of scaling with
applications to the 1968 and 1972 presidential elections. *American
Political Science Review*, 71(1), 111-130.
[doi:10.2307/1956957](https://doi.org/10.2307/1956957)

Poole, K. T. (1998). Recovering a basic space from a set of issue
scales. *American Journal of Political Science*, 42(3), 954-993.
[doi:10.2307/2991737](https://doi.org/10.2307/2991737)

## See also

Other scaling:
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
[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md),
[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)

## Examples

``` r
# 40 simulated respondents place five stimuli on a 1-7 scale
set.seed(1)
truth <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
Z <- t(vapply(1:40, function(i) {
  4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * truth + rnorm(5, 0, 0.3)
}, numeric(5)))
round(amscal(Z)$zhat, 2)
#> [1] -1.28 -0.59 -0.05  0.67  1.25
round((truth - mean(truth)) / sd(truth), 2)
#> [1] -1.30 -0.55 -0.08  0.68  1.24
```
