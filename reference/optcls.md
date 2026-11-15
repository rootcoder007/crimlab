# Optimal Classification scaling

Nonparametric ideal points from binary votes in the spirit of Poole
(2000): each vote's cutting direction is the unit-length difference of
the mean positions of its yea and nay voters, each legislator moves to
the mean of the directions of the sides taken, for `maxit` sweeps from
each of `nstart` random starts; the start with the fewest classification
errors is kept. A vote is split at the point of its direction that
misclassifies fewest of its voters, either side voting yea, so no vote
has more errors than voting with the majority.

## Usage

``` r
optcls(votes, ndim = 1, maxit = 500, nstart = 10, seed = 42)
```

## Arguments

- votes:

  A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.

- ndim:

  The number of dimensions.

- maxit:

  The sweeps per start.

- nstart:

  The random starts.

- seed:

  The seed of the starts, as
  [`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

## Value

A list of class `optcls`: `ideal`, `normals` (votes by dimensions),
`pre` (the proportional reduction in error), `errors`, `null` (the
errors of voting with the majority on every vote) and `ndim`.

## References

Poole, K. T. (2000). Nonparametric unfolding of binary choice data.
*Political Analysis*, 8(3), 211-237.
[doi:10.1093/oxfordjournals.pan.a029814](https://doi.org/10.1093/oxfordjournals.pan.a029814)

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
# 20 simulated legislators, 30 votes cut along one dimension
set.seed(1)
x <- sort(rnorm(20))
v <- outer(x, rnorm(30), ">") * 1
optcls(v, maxit = 20, nstart = 2)$pre
#> [1] 1
```
