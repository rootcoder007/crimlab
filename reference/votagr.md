# Vote agreement between legislators

The share of the roll calls on which each pair of legislators (judges,
council members) voted the same way, over the roll calls both took part
in (Armstrong et al. 2021, section 3.2.2): the similarity matrix that
metric scaling of a chamber starts from.

## Usage

``` r
votagr(x)
```

## Arguments

- x:

  A legislator-by-vote matrix of votes (any coding; equal values agree),
  `NA` for an absence. A vector is one vote.

## Value

A list of class `votagr`: `agree` (legislators by legislators, 1 on the
diagonal, `NA` for a pair with no vote in common), `mean` (the mean
agreement over the pairs), `n` and `m` (legislators and votes).

## References

Armstrong, D. A., Bakker, R., Carroll, R., Hare, C., Poole, K. T. and
Rosenthal, H. (2021). *Analyzing Spatial Models of Choice and Judgment*,
2nd ed. Chapman & Hall/CRC.

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
[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)

## Examples

``` r
# five judges on eight simulated cases
set.seed(1)
v <- matrix(rbinom(40, 1, 0.5), 5, 8)
v[2, 3] <- NA
votagr(v)$agree
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.5714286 0.5000000 0.6250000 0.2500000
#> [2,] 0.5714286 1.0000000 0.4285714 0.4285714 0.7142857
#> [3,] 0.5000000 0.4285714 1.0000000 0.3750000 0.5000000
#> [4,] 0.6250000 0.4285714 0.3750000 1.0000000 0.1250000
#> [5,] 0.2500000 0.7142857 0.5000000 0.1250000 1.0000000
```
