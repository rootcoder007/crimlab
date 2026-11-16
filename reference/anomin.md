# Alpha-NOMINATE

Legislator \\i\\ votes yea on roll call \\j\\ with probability
\\\Phi(u\_{ij})\\, \\u\_{ij} = Q\_{ij} + \alpha (G\_{ij} - Q\_{ij})\\,
where \\Q\\ is the quadratic and \\G\\ the Gaussian utility difference
of the distances to the yea and nay outcomes, with weight 0.5 (Carroll
et al. 2013): \\\alpha = 1\\ is Gaussian utility, \\\alpha = 0\\
quadratic. \\\beta\\ has a flat prior on \\(0, \infty)\\, \\\alpha\\ on
\\\[0, 1\]\\; coordinates carry the hierarchical penalty \\-x'Sx/2\\
with \\S\\ redrawn each sweep from an inverse Wishart. Every coordinate,
\\\beta\\ and \\\alpha\\ are slice sampled (Neal 2003: width 8, at most
3 steps out); starts are uniform on \\(-1, 1)\\, \\\beta = 10\\,
\\\alpha = 0.7\\. Each draw is then centred on its legislators and
rotated onto the posterior mean, the `polarity` legislator positive.
Roll calls whose minority side is below `lop` and legislators with fewer
than `minvotes` votes are dropped first. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Usage

``` r
anomin(
  votes,
  ndim = 1,
  nsamp = 500,
  burn = 100,
  seed = 42,
  thin = 1,
  lop = 0.025,
  minvotes = 20,
  polarity = 1,
  fixalpha = FALSE
)
```

## Arguments

- votes:

  A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.

- ndim:

  The number of dimensions; `ndim` times the points scaled must be at
  most 67108864 (2^26).

- nsamp:

  The draws kept; `nsamp` times the values kept per draw must be at most
  67108864 (2^26).

- burn:

  The sweeps discarded first; `burn` plus `nsamp` times `thin` must be
  at most 2147483647.

- seed:

  The seed.

- thin:

  Keep every `thin`-th sweep after `burn`.

- lop:

  The smallest minority share of a kept roll call.

- minvotes:

  The fewest kept roll calls of a kept legislator.

- polarity:

  The row of the legislator placed on the positive side of each
  dimension (recycled).

- fixalpha:

  Whether to fix \\\alpha = 1\\.

## Value

A list of class `anomin`: `ideal` and `sd` (kept legislators by
dimensions), `yea` and `nay` (outcome locations), `alpha`, `interval`
(2.5% and 97.5% of \\\alpha\\), `beta`, `draws` (`X`, `Y`, `N` arrays,
draws first, and `beta`, `alpha`), `legis` and `votes` (the kept rows
and columns) and `ndim`.

## References

Carroll, R., Lewis, J. B., Lo, J., Poole, K. T. and Rosenthal, H.
(2013). The structure of utility in spatial models of voting. *American
Journal of Political Science*, 57(4), 1008-1028.
[doi:10.1111/ajps.12029](https://doi.org/10.1111/ajps.12029)

Neal, R. M. (2003). Slice sampling. *Annals of Statistics*, 31(3),
705-767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
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
# 20 simulated legislators on a line, 30 roll calls, some votes flipped
set.seed(1)
x <- seq(-1, 1, length.out = 20)
V <- outer(x, runif(30, -0.7, 0.7), ">") * 1
flip <- sample(length(V), 30)
V[flip] <- 1 - V[flip]
f <- anomin(V, nsamp = 20, burn = 10, minvotes = 10, polarity = 20)
f$alpha
#> [1] 0.9071232
```
