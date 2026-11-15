# NOMINATE utilities and vote probabilities

The Gaussian utility of Poole and Rosenthal (1985): legislator \\i\\
gains \\U = \beta \exp(-\frac{1}{2} \sum_k w_k^2 (x\_{ik} -
z\_{jk})^2)\\ from the outcome \\z_j\\; the probability of a yea is the
logistic function of the yea utility minus the nay utility.

## Usage

``` r
nomutl(x, zy, zn, beta = 15, w = NULL)
```

## Arguments

- x:

  Legislators' ideal points (legislators by dimensions).

- zy:

  The yea outcomes (votes by dimensions).

- zn:

  The nay outcomes (votes by dimensions).

- beta:

  The signal-to-noise ratio.

- w:

  `NULL` for weights of 1, or one weight per dimension.

## Value

A list of class `nomutl`: `uyea`, `unay`, `udiff` and `prob`
(legislators by votes).

## References

Poole, K. T. and Rosenthal, H. (1985). A spatial model for legislative
roll call analysis. *American Journal of Political Science*, 29(2),
357-384. [doi:10.2307/2111172](https://doi.org/10.2307/2111172)

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
nomutl(matrix(rnorm(8), 4), matrix(rnorm(6), 3), matrix(rnorm(6), 3))$prob
#>             [,1]       [,2]         [,3]
#> [1,] 0.989937668 0.10357690 6.370307e-06
#> [2,] 0.990493247 0.99990393 4.484893e-02
#> [3,] 0.971693451 0.03457389 2.401618e-05
#> [4,] 0.008546933 0.05135459 2.061332e-02
```
