# Bayesian Aldrich-McKelvey scaling

Respondent \\i\\ places stimulus \\j\\ at \\z\_{ij} \sim N(a_i + b_i
\hat z_j, 1 / (\tau_i \tau_j))\\: a respondent's shift and stretch of
the common positions \\\hat z\\, with noise whose precision factors into
a respondent and a stimulus part (Hare et al. 2015). Priors: \\a_i, b_i
\sim U(-100, 100)\\, \\\tau_j \sim G(0.1, 0.1)\\, \\\tau_i \sim G(g_a,
g_b)\\ with \\g_a, g_b \sim G(0.1, 0.1)\\, and \\\hat z\\ the
standardised \\z^\* \sim N(0, 1)\\, the `polarity` stimulus restricted
to the left. Gibbs sampling: \\a_i\\ and \\b_i\\ from truncated normals,
the precisions and \\g_b\\ from conjugate gammas, \\g_a\\ and each
\\z^\*\_j\\ by slice sampling (Neal 2003). Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it; rows with
no placement are dropped.

## Usage

``` r
bamscl(Z, nsamp = 1000, burn = 200, polarity = 1, seed = 42)
```

## Arguments

- Z:

  A respondent-by-stimulus matrix of placements, `NA` for missing.

- nsamp:

  The draws kept; `nsamp` times the values kept per draw must be at most
  67108864 (2^26).

- burn:

  The sweeps discarded first; `burn` plus `nsamp` must be at most
  2147483647.

- polarity:

  The column of the stimulus placed on the left.

- seed:

  The seed.

## Value

A list of class `bamscl`: `zeta`, `sd` and `interval` (2.5% and 97.5%
rows) of the stimulus positions, `a` and `b` (posterior means), `tau`
(stimulus precisions), `draws` (draws by stimuli) and `nsamp`.

## References

Hare, C., Armstrong, D. A., Bakker, R., Carroll, R. and Poole, K. T.
(2015). Using Bayesian Aldrich-McKelvey scaling to study citizens'
ideological preferences and perceptions. *American Journal of Political
Science*, 59(3), 759-774.
[doi:10.1111/ajps.12151](https://doi.org/10.1111/ajps.12151)

Neal, R. M. (2003). Slice sampling. *Annals of Statistics*, 31(3),
705-767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

## See also

Other scaling:
[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md),
[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md),
[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md),
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
# 40 simulated respondents place four stimuli
set.seed(1)
zt <- c(-1, -0.3, 0.4, 1.2)
Z <- rnorm(40, 0, 0.5) + outer(runif(40, 0.5, 1.5), zt) +
  matrix(rnorm(160, 0, 0.3), 40, 4)
round(bamscl(Z, nsamp = 100, burn = 50)$zeta, 2)
#> [1] -1.14 -0.42  0.40  1.16
```
