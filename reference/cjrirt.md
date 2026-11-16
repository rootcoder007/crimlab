# Clinton-Jackman-Rivers Bayesian IRT

The probit ideal-point model \\P(y\_{ij} = 1) = \Phi(\beta_j' x_i -
\alpha_j)\\ with priors \\x_i \sim N(0, I)\\ and \\(\alpha_j, \beta_j)
\sim N(0, v I)\\ (Clinton, Jackman and Rivers 2004), by the Gibbs
sampler of Albert and Chib (1993): latent utilities from truncated
normals, then every \\(\alpha_j, \beta_j)\\ and every \\x_i\\ from its
normal full conditional. Every draw is an inverse normal of a uniform of
Philox4x32-10 (Salmon et al. 2011), sweep \\t\\ on stream \\t\\ of key
`seed`, so a chain is the same on every platform. The default start is
the posterior mode of
[`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md).

## Usage

``` r
cjrirt(
  votes,
  ndim = 1,
  nsamp = 1000,
  burn = 200,
  bvar = 25,
  start = NULL,
  seed = 0,
  chains = TRUE
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

  The sweeps discarded first; `burn` plus `nsamp` must be at most
  2147483647.

- bvar:

  The prior variance \\v\\ of the vote parameters.

- start:

  `NULL`, or starting ideal points.

- seed:

  The Philox key, a whole number from 0.

- chains:

  Whether to return the draws as well.

## Value

A list of class `cjrirt`: `ideal`, `sd`, `alpha`, `beta`, `nsamp` and,
with `chains`, `xchain` (draws by legislators by dimensions), `achain`
and `bchain`.

## References

Clinton, J., Jackman, S. and Rivers, D. (2004). The statistical analysis
of roll call data. *American Political Science Review*, 98(2), 355-370.
[doi:10.1017/S0003055404001194](https://doi.org/10.1017/S0003055404001194)

Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669-679.
[doi:10.1080/01621459.1993.10476321](https://doi.org/10.1080/01621459.1993.10476321)

Salmon, J. K., Moraes, M. A., Dror, R. O. and Shaw, D. E. (2011).
Parallel random numbers: as easy as 1, 2, 3. *Proceedings of the
International Conference for High Performance Computing, Networking,
Storage and Analysis* (SC11).
[doi:10.1145/2063384.2063405](https://doi.org/10.1145/2063384.2063405)

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md),
[`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md),
[`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md),
[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md),
[`bunfld()`](https://rootcoder007.github.io/crimlab/reference/bunfld.md),
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
votes <- rbind(c(1, 1, 0, 0), c(1, 1, 1, 0), c(0, 1, 1, 1),
               c(0, 0, 1, 1), c(1, 0, 0, 0))
round(cjrirt(votes, nsamp = 50, burn = 10)$ideal, 3)
#>        [,1]
#> [1,] -1.209
#> [2,] -0.508
#> [3,]  0.358
#> [4,]  0.859
#> [5,] -0.904
```
