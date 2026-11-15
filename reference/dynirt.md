# Dynamic IRT (ideal points that drift)

The one-dimensional dynamic item response model of Martin and Quinn
(2002), made for courts: a vote is cast when \\z\_{jk} = -\alpha_k +
\beta_k \theta\_{j,t(k)} + \varepsilon\_{jk} \> 0\\; the ideal points
follow a random walk, \\\theta\_{j,0} \sim N(e_0, E_0)\\,
\\\theta\_{j,t} \sim N(\theta\_{j,t-1}, \tau^2_j)\\; \\\alpha_k \sim
N(a_0, 1/A_0)\\, \\\beta_k \sim N(b_0, 1/B_0)\\. Each sweep draws the
latent utilities, every case's \\(\alpha, \beta)\\, every justice's path
by forward filtering and backward sampling and, when `c0` and `d0` are
positive, the evolution variances from \\IG((c_0 + T)/2, (d_0 +
SS)/2)\\. A draw is reflected whenever the `anchor` member's mean ideal
point is negative. The chain starts from the first principal component
of the votes. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
dynirt(
  votes,
  period,
  nsamp = 500,
  burn = 100,
  seed = 42,
  thin = 1,
  tau2 = 1,
  c0 = -1,
  d0 = -1,
  e0 = 0,
  E0 = 1,
  a0 = 0,
  A0 = 0.1,
  b0 = 0,
  B0 = 0.1,
  anchor = NULL
)
```

## Arguments

- votes:

  A member-by-vote matrix of 0/1 votes (-1/1 is read as 0/1), `NA`
  absent.

- period:

  The period (term) of each vote, one per column.

- nsamp:

  The draws kept.

- burn:

  The sweeps discarded first.

- seed:

  The seed.

- thin:

  Keep every `thin`-th sweep.

- tau2:

  The evolution variance (its start when estimated).

- c0, d0:

  The inverse-gamma prior of \\\tau^2\\, held at `tau2` unless both are
  positive.

- e0, E0:

  The prior mean and variance of the first ideal point.

- a0, A0:

  The prior mean and precision of \\\alpha\\.

- b0, B0:

  The prior mean and precision of \\\beta\\.

- anchor:

  `NULL` (the largest first principal component score) or the row kept
  on the positive side.

## Value

A list of class `dynirt`: `theta` and `sd` (members by periods),
`alpha`, `beta`, `tau2`, `periods`, `nsamp` and `anchor`.

## References

Martin, A. D. and Quinn, K. M. (2002). Dynamic ideal point estimation
via Markov chain Monte Carlo for the U.S. Supreme Court, 1953-1999.
*Political Analysis*, 10(2), 134-153.
[doi:10.1093/pan/10.2.134](https://doi.org/10.1093/pan/10.2.134)

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
# 12 simulated justices over three terms of 15 cases
set.seed(7)
period <- rep(1:3, each = 15)
walk <- function(s) s + cumsum(c(0, rnorm(2, 0, 0.3)))
truth <- t(sapply(rnorm(12), walk))
p <- pnorm(-rep(rnorm(45, 0, 0.5), each = 12) +
  rep(rnorm(45, 1.5, 0.5), each = 12) * truth[, period])
V <- matrix(rbinom(length(p), 1, p), 12)
f <- dynirt(V, period, nsamp = 100, burn = 50)
cor(as.numeric(f$theta), as.numeric(truth))
#> [1] 0.9013526
```
