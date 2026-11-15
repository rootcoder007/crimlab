# EM estimation of binary IRT ideal points

The EM algorithm of Imai, Lo and Olmsted (2016) for the probit model
\\y^\*\_{ij} = \alpha_j + \beta_j' x_i + e\_{ij}\\, yea when
\\y^\*\_{ij} \> 0\\, with priors \\x_i \sim N(\mu_x, s^2_x I)\\ and
\\(\alpha_j, \beta_j) \sim N(\mu_b, s^2_b I)\\. The E-step replaces
\\y^\*\\ by its truncated-normal mean (untruncated when absent); the
M-step is the pair of ridge regressions; the fixed point is the
posterior mode. The default start is \\\alpha = \beta = 0\\ and \\x\\
from 50 steps of orthogonal iteration on the cross-product of the
column-centred +1/-1 votes.

## Usage

``` r
irtem(
  votes,
  ndim = 1,
  maxit = 500,
  tol = 1e-06,
  xprior = c(0, 1),
  bprior = c(0, 25),
  start = NULL,
  conv = "cor"
)
```

## Arguments

- votes:

  A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.

- ndim:

  The number of dimensions.

- maxit:

  The most iterations.

- tol:

  The convergence threshold.

- xprior:

  The prior mean and variance of the ideal points.

- bprior:

  The prior mean and variance of the intercepts and slopes.

- start:

  `NULL`, or a list of the starting intercepts, slopes and ideal points.

- conv:

  `"cor"` (1 minus the smallest correlation of old and new values,
  column by column) or `"abs"` (the largest absolute change).

## Value

A list of class `irtem`: `ideal`, `disc`, `diff` (the intercepts),
`loglik`, `iter` and `conv`.

## References

Imai, K., Lo, J. and Olmsted, J. (2016). Fast estimation of ideal points
with massive data. *American Political Science Review*, 110(4), 631-656.
[doi:10.1017/S000305541600037X](https://doi.org/10.1017/S000305541600037X)

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
x <- rnorm(20)
v <- matrix(rbinom(600, 1, pnorm(outer(x, rnorm(30)))), 20)
cor(irtem(v)$ideal[, 1], x)
#> [1] 0.8854776
```
