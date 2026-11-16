# Parametric bootstrap of ideal points

Standard errors of ideal points by the parametric bootstrap of Lewis and
Poole (2004): votes are simulated from the fitted logistic probabilities
(\\\beta = 15\\ times the projection minus the cutpoint), absences kept,
and each legislator re-estimated by 20 gradient steps from the fitted
point plus normal noise of standard deviation 0.1. Draws come from
`seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
nombt(votes, x, normals, cuts, nboot = 100, seed = 42)
```

## Arguments

- votes:

  A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.

- x:

  Fitted ideal points.

- normals:

  Fitted normal vectors (votes by dimensions).

- cuts:

  Fitted cutpoints, one per vote.

- nboot:

  The bootstrap replicates; `nboot` times the ideal points must be at
  most 67108864 (2^26).

- seed:

  The seed.

## Value

A list of class `nombt`: `se` and `mean` (legislators by dimensions) and
`nboot`.

## References

Lewis, J. B. and Poole, K. T. (2004). Measuring bias and uncertainty in
ideal point estimates via the parametric bootstrap. *Political
Analysis*, 12(2), 105-127.

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
v <- matrix(rbinom(40, 1, 0.5), 5, 8)
f <- dwnom(v, maxit = 5)
nombt(v, f$ideal, f$normals, f$cuts, nboot = 5)$se
#>            [,1]       [,2]
#> [1,] 0.03927144 0.08696910
#> [2,] 0.07038480 0.10135607
#> [3,] 0.08214763 0.11470586
#> [4,] 0.06887660 0.13162676
#> [5,] 0.13492223 0.07540698
```
