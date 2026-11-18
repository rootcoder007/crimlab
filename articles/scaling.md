# Scaling and ideal points

Agreement scores and the scaling methods that place voters, judges or
items on a line or in a plane.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Vote agreement between legislators

[`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md):
The share of the roll calls on which each pair of legislators (judges,
council members) voted the same way, over the roll calls both took part
in (Armstrong et al. 2021, section 3.2.2): the similarity matrix that
metric scaling of a chamber starts from.

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

## Aldrich-McKelvey scaling

[`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md):
Recovers the positions of stimuli (parties, candidates, courts) on one
scale from respondents’ placements of them, each respondent allowed an
own shift and stretch of the scale: z\_{ij} = a_i + b_i z_j +
epsilon\_{ij}. The positions are the least-squares solution of Aldrich
and McKelvey (1977): with X_i = \[1, z_i\] and P_i the projection onto
its columns, z is the eigenvector of sum_i (I - P_i) with the smallest
eigenvalue orthogonal to the constant. Respondents with a missing
placement, or who place every stimulus at the same point, are left out
of that sum. The positions are standardised to mean 0 and standard
deviation 1, the first stimulus on the left. The respondent intercepts
and slopes are each respondent’s least-squares regression of the
placements on those positions (defined with two placements or more; a
slope of zero is set to 1e-10).

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

## Basic Space (blackbox) scaling

[`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md):
Recovers respondents’ positions from a respondent-by-issue matrix of
responses by the decomposition of Poole (1998), X_0 = Psi W’ + J_n c’ +
E_0 over the observed cells. Each dimension starts from the
sign-corrected respondent means (column polarities from the pairwise
correlations) and is refined by four alternating least-squares sweeps on
the residuals of the earlier dimensions; up to five joint sweeps follow,
stopping once the two half-step error sums of squares differ by less
than 0.01; the fitted Psi W’ is then re-expressed by its singular value
decomposition, Psi = U D^{1/2}, W = V D^{1/2}. Respondents with fewer
than `minscale` responses are not scaled. The sign of each dimension is
arbitrary.

``` r

# 30 simulated respondents on eight issues, ten answers missing
set.seed(1)
X <- matrix(rnorm(60), 30, 2) %*% matrix(rnorm(16), 2, 8) +
  matrix(rnorm(240, 0, 0.3), 30, 8)
X[sample(240, 10)] <- NA
bbscal(X, ndim = 2)$fit
#> [1] 0.9787813
```

## Optimal Classification scaling

[`optcls()`](https://rootcoder007.github.io/crimlab/reference/optcls.md):
Nonparametric ideal points from binary votes in the spirit of Poole
(2000): each vote’s cutting direction is the unit-length difference of
the mean positions of its yea and nay voters, each legislator moves to
the mean of the directions of the sides taken, for `maxit` sweeps from
each of `nstart` random starts; the start with the fewest classification
errors is kept. A vote is split at the point of its direction that
misclassifies fewest of its voters, either side voting yea, so no vote
has more errors than voting with the majority.

``` r

# 20 simulated legislators, 30 votes cut along one dimension
set.seed(1)
x <- sort(rnorm(20))
v <- outer(x, rnorm(30), ">") * 1
optcls(v, maxit = 20, nstart = 2)$pre
#> [1] 1
```

## Ordered Optimal Classification

[`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md):
Ordinal responses (a survey’s agree-disagree items, a court’s graded
outcomes) classified by item directions and cutpoints in the spirit of
Hare, Liu and Lupton (2018): respondents and item directions start
normal (seed 42, as [`set.seed()`](https://rdrr.io/r/base/Random.html)
takes it), the cutpoints of each item at equal steps on \[-1, 1\]; each
sweep moves every cutpoint to the median of the midpoints between
neighbouring respondents (in projection order) on either side of it,
until the share classified correctly stops changing.

``` r

set.seed(1)
Y <- matrix(sample(1:4, 60, replace = TRUE), 15, 4)
ordocl(Y, ndim = 1, maxit = 20)$correct
#> [1] 0.2
```

## Double-centre a distance matrix

[`dblctr()`](https://rootcoder007.github.io/crimlab/reference/dblctr.md):
B = - H D^{(2)} H with H = I - n^{-1} 1 1’ and D^{(2)} the squared
distances: the inner-product matrix of classical scaling (Torgerson
1952).

``` r

D <- as.matrix(dist(cbind(c(0, 1, 3, 4), c(0, 1, 0, 2))))
dblctr(D)$B
#>         [,1]    [,2]    [,3]    [,4]
#> [1,]  4.5625  1.8125 -1.4375 -4.9375
#> [2,]  1.8125  1.0625 -1.1875 -1.6875
#> [3,] -1.4375 -1.1875  1.5625  1.0625
#> [4,] -4.9375 -1.6875  1.0625  5.5625
```

## Classical (metric) multidimensional scaling

[`cmdsf()`](https://rootcoder007.github.io/crimlab/reference/cmdsf.md):
Torgerson scaling: the leading eigenvectors of the double-centred
squared distances
([`dblctr()`](https://rootcoder007.github.io/crimlab/reference/dblctr.md)),
each scaled by the square root of its eigenvalue (negative ones set to
0). The fit is Mardia’s criterion, the kept positive eigenvalues over
the sum of all absolute ones; the stress is Kruskal’s stress-1 against
the positive distances. The sign of each dimension is arbitrary.

``` r

set.seed(1)
D <- as.matrix(dist(matrix(rnorm(40), 10)))
cmdsf(D)$fit
#> [1] 0.7816266
```

## SMACOF stress minimisation

[`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md):
Metric scaling by majorisation (de Leeuw 1977): each step is the Guttman
transform X arrow V^{+} B(X) X with V = (W 1) - W, until the raw stress
sum w\_{ij} (d\_{ij} - delta\_{ij})^2 changes by less than `tol`.

``` r

set.seed(1)
D <- as.matrix(dist(matrix(rnorm(40), 10)))
smacf(D)$stress
#> [1] 6.626513
```

## Nonmetric multidimensional scaling

[`nmdsf()`](https://rootcoder007.github.io/crimlab/reference/nmdsf.md):
Kruskal’s nonmetric scaling: at each step the distances of the current
configuration are fitted by a monotone (isotonic, pool-adjacent-
violators) regression on the order of the dissimilarities, and one
SMACOF step
([`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md))
is taken towards the fitted disparities, until no coordinate moves by
`tol`. The stress is Kruskal’s stress-1 of the last step. The start is
normal from seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

``` r

set.seed(1)
D <- as.matrix(dist(matrix(rnorm(40), 10)))
nmdsf(D, maxit = 50)$stress
#> [1] 0.1259758
```

## Fit of a scaling by its eigenvalues

[`mdsfit()`](https://rootcoder007.github.io/crimlab/reference/mdsfit.md):
Mardia’s criterion: the share of each dimension, its eigenvalue (0 when
negative) over the sum of the absolute values of all eigenvalues, and
the cumulative shares.

``` r

mdsfit(c(4, 2, 1, 0.5))
#> mdsfit:
#>   fit    0.53330 0.26670 0.13330 0.06667
#>   cumfit 0.5333 0.8000 0.9333 1.0000
#>   ev     4.0 2.0 1.0 0.5
```

## INDSCAL: scaling with individual differences

[`indscl()`](https://rootcoder007.github.io/crimlab/reference/indscl.md):
The weighted Euclidean model of Carroll and Chang (1970): one stimulus
space shared by every individual, each individual weighting its
dimensions. The configuration starts from classical scaling of the mean
dissimilarities; then the weights of each individual and dimension (at
least 0.01) and the stimulus coordinates are updated in turn until the
configuration changes by a relative `tol`. The stress is the sum of
squared differences between each individual’s dissimilarities and the
distances in its weighted space.

``` r

set.seed(1)
D1 <- as.matrix(dist(matrix(rnorm(20), 5)))
D2 <- as.matrix(dist(matrix(rnorm(20), 5)))
indscl(list(D1, D2), maxit = 30)$weights
#>      [,1]    [,2]
#> [1,]   NA 1979597
#> [2,]   NA 1929755
```

## Stress of an unfolding

[`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md):
The weighted sum of squared differences between the distances from
respondents to stimuli in a configuration and the observed
dissimilarities, over the observed cells.

``` r

set.seed(1)
unfstr(matrix(rnorm(6), 3), matrix(rnorm(8), 4), matrix(runif(12), 3))
#> unfstr:
#>   stress 32.263
```

## Least-squares metric unfolding

[`mlsmu()`](https://rootcoder007.github.io/crimlab/reference/mlsmu.md):
Unfolding of a respondent-by-stimulus matrix in the alternating
least-squares form of MLSMU6 (Poole 1984): respondents and stimuli take
gradient steps in turn against the row-centred dissimilarities, each
step scaled by the largest eigenvalue of the gradient’s cross-product
and both clouds centred, until the stress
([`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md))
changes by a relative `tol`; of `nstart` normal starts (seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it) the lowest
stress is kept.

``` r

set.seed(1)
D <- matrix(runif(120), 20, 6)
mlsmu(D, nstart = 1, maxit = 50)$stress
#> [1] 111.1204
```

## SMACOF unfolding

[`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md):
Unfolding by majorisation: respondents and stimuli are embedded in one
space by SMACOF steps
([`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md))
on the joint configuration with weight 1 between a respondent and a
stimulus and 0 within each set (Borg and Groenen 2005, chapter 14). The
start is normal from seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

``` r

set.seed(1)
smunf(matrix(runif(12), 3), maxit = 20)$stress
#> [1] 3.242821
```

## Ideal points from an unfolding

[`ipsrec()`](https://rootcoder007.github.io/crimlab/reference/ipsrec.md):
In an unfolding the respondents’ coordinates are their ideal points
(Armstrong et al. 2021, section 4.5): the respondent coordinates as a
matrix.

``` r

set.seed(1)
u <- smunf(matrix(runif(12), 3), maxit = 20)
ipsrec(u$xr, u$xs)$ideal
#>            [,1]       [,2]
#> [1,] -0.4765297 -0.2645617
#> [2,] -0.3149654 -0.1746451
#> [3,] -0.4975857 -0.2764589
```

## Procrustes rotation

[`procrs()`](https://rootcoder007.github.io/crimlab/reference/procrs.md):
The orthogonal rotation (a proper rotation: no reflection) that best
carries the centred configuration `X` onto the centred `target` in least
squares, the target’s centre added back (Gower and Dijksterhuis 2004).

``` r

set.seed(1)
A <- matrix(rnorm(20), 10)
procrs(A, A %*% matrix(c(0, 1, -1, 0), 2))$rotation
#>      [,1]          [,2]
#> [1,]    0 -1.000000e+00
#> [2,]    1 -1.110223e-16
```

## NOMINATE utilities and vote probabilities

[`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md):
The Gaussian utility of Poole and Rosenthal (1985): legislator i gains U
= beta (- sum_k w_k^2 (x\_{ik} - z\_{jk})^2) from the outcome z_j; the
probability of a yea is the logistic function of the yea utility minus
the nay utility.

``` r

set.seed(1)
nomutl(matrix(rnorm(8), 4), matrix(rnorm(6), 3), matrix(rnorm(6), 3))$prob
#>             [,1]       [,2]         [,3]
#> [1,] 0.989937668 0.10357690 6.370307e-06
#> [2,] 0.990493247 0.99990393 4.484893e-02
#> [3,] 0.971693451 0.03457389 2.401618e-05
#> [4,] 0.008546933 0.05135459 2.061332e-02
```

## One NOMINATE vote probability

[`nompr()`](https://rootcoder007.github.io/crimlab/reference/nompr.md):
The probability of a yea of
[`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md)
for one legislator and one vote.

``` r

nompr(0.1, 0.3, -0.3)$prob
#> [1] 0.7018734
```

## NOMINATE log-likelihood and classification

[`nomll()`](https://rootcoder007.github.io/crimlab/reference/nomll.md):
The log-likelihood of the observed votes under the probabilities of
[`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md)
(clamped to \[10^{-10}, 1 - 10^{-10}\]) and `gmp`, the share of votes
the model classifies correctly.

``` r

set.seed(1)
v <- matrix(rbinom(20, 1, 0.5), 4, 5)
nomll(v, matrix(rnorm(4)), matrix(rnorm(5)), matrix(rnorm(5)))
#> nomll:
#>   loglik  -48.178
#>   gmp     0.6
#>   correct 12
#>   total   20
```

## DW-NOMINATE-style ideal points

[`dwnom()`](https://rootcoder007.github.io/crimlab/reference/dwnom.md):
A deterministic fixed-point fit in the Gaussian-utility family of Poole
and Rosenthal (1997): from normal starts (seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it), each
vote’s normal vector and midpoint follow the mean positions of its two
sides, and each legislator moves to the mean of the outcome locations
(midpoint plus or minus half the normal) of the sides taken, the cloud
shrunk into the unit ball, for `maxit` sweeps. The log-likelihood and
the share classified correctly use the probit of the Gaussian utility
difference with beta = 15 and equal dimension weights.

``` r

set.seed(1)
v <- matrix(rbinom(600, 1, 0.5), 20, 30)
dwnom(v, maxit = 20)$gmp
#> [1] 0.7033333
```

## Parametric bootstrap of ideal points

[`nombt()`](https://rootcoder007.github.io/crimlab/reference/nombt.md):
Standard errors of ideal points by the parametric bootstrap of Lewis and
Poole (2004): votes are simulated from the fitted logistic probabilities
(beta = 15 times the projection minus the cutpoint), absences kept, and
each legislator re-estimated by 20 gradient steps from the fitted point
plus normal noise of standard deviation 0.1. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

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

## Binary IRT likelihood

[`irtlik()`](https://rootcoder007.github.io/crimlab/reference/irtlik.md):
The log-likelihood of the votes under the probit item response model of
Clinton, Jackman and Rivers (2004), P(y\_{ij} = 1) = Phi(beta_j’ x_i -
alpha_j) (clamped to \[10^{-10}, 1 - 10^{-10}\]), with the share of
votes classified correctly.

``` r

set.seed(1)
v <- matrix(rbinom(20, 1, 0.5), 4, 5)
irtlik(v, matrix(rnorm(4)), rep(0, 5), matrix(rnorm(5)))$loglik
#> [1] -22.8751
```

## EM estimation of binary IRT ideal points

[`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md):
The EM algorithm of Imai, Lo and Olmsted (2016) for the probit model
y^**{ij} = alpha_j + beta_j’ x_i + e*{ij}, yea when y^*\_{ij} \> 0, with
priors x_i N(mu_x, s^2_x I) and (alpha_j, beta_j) N(mu_b, s^2_b I). The
E-step replaces y^\* by its truncated-normal mean (untruncated when
absent); the M-step is the pair of ridge regressions; the fixed point is
the posterior mode. The default start is alpha = beta = 0 and x from 50
steps of orthogonal iteration on the cross-product of the column-centred
+1/-1 votes.

``` r

set.seed(1)
x <- rnorm(20)
v <- matrix(rbinom(600, 1, pnorm(outer(x, rnorm(30)))), 20)
cor(irtem(v)$ideal[, 1], x)
#> [1] 0.8854776
```

## Posterior summaries of ideal-point draws

[`irtpst()`](https://rootcoder007.github.io/crimlab/reference/irtpst.md):
Means, standard deviations and central 95% intervals (R’s default
quantile, type 7) of a chain of ideal points, each draw first
standardised dimension by dimension when `std` is `TRUE` (location and
scale are not identified in IRT models; Jackman 2009).

``` r

set.seed(1)
irtpst(array(rnorm(100 * 5 * 2), c(100, 5, 2)))$mean
#>             [,1]         [,2]
#> [1,]  0.06747061 -0.008652293
#> [2,] -0.06792212 -0.144214340
#> [3,]  0.01295946  0.064296854
#> [4,]  0.07157441  0.079645425
#> [5,] -0.08408237  0.008924354
```

## Clinton-Jackman-Rivers Bayesian IRT

[`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md):
The probit ideal-point model P(y\_{ij} = 1) = Phi(beta_j’ x_i - alpha_j)
with priors x_i N(0, I) and (alpha_j, beta_j) N(0, v I) (Clinton,
Jackman and Rivers 2004), by the Gibbs sampler of Albert and Chib
(1993): latent utilities from truncated normals, then every (alpha_j,
beta_j) and every x_i from its normal full conditional. Every draw is an
inverse normal of a uniform of Philox4x32-10 (Salmon et al. 2011), sweep
t on stream t of key `seed`, so a chain is the same on every platform.
The default start is the posterior mode of
[`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md).

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

## Bayesian Aldrich-McKelvey scaling

[`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md):
Respondent i places stimulus j at z\_{ij} N(a_i + b_i z_j, 1 / (tau_i
tau_j)): a respondent’s shift and stretch of the common positions z,
with noise whose precision factors into a respondent and a stimulus part
(Hare et al. 2015). Priors: a_i, b_i U(-100, 100), tau_j G(0.1, 0.1),
tau_i G(g_a, g_b) with g_a, g_b G(0.1, 0.1), and z the standardised z^\*
N(0, 1), the `polarity` stimulus restricted to the left. Gibbs sampling:
a_i and b_i from truncated normals, the precisions and g_b from
conjugate gammas, g_a and each z^\*\_j by slice sampling (Neal 2003).
Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it; rows with
no placement are dropped.

``` r

# 40 simulated respondents place four stimuli
set.seed(1)
zt <- c(-1, -0.3, 0.4, 1.2)
Z <- rnorm(40, 0, 0.5) + outer(runif(40, 0.5, 1.5), zt) +
  matrix(rnorm(160, 0, 0.3), 40, 4)
round(bamscl(Z, nsamp = 100, burn = 50)$zeta, 2)
#> [1] -1.14 -0.42  0.40  1.16
```

## Bayesian metric multidimensional scaling

[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md):
Each observed dissimilarity is lognormal about the distance between the
two objects’ coordinates, delta\_{ij} N(d\_{ij}, 1 / tau), with x\_{ik}
N(0, 10^2) and tau U(0, 10) (Bakker and Poole 2013). Coordinates are
slice sampled one at a time (Neal 2003), tau drawn from its truncated
gamma; the chain starts from classical scaling plus small normal noise.
The likelihood is unchanged by translation, rotation and reflection, so
every draw is aligned (translation and rotation) onto the posterior
mean. Zero, negative or missing dissimilarities are left out. Draws come
from `seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html) takes
it.

``` r

set.seed(1)
D <- as.matrix(dist(matrix(rnorm(16), 8)))
f <- bmdsf(D, nsamp = 50, burn = 50)
cor(f$dist[upper.tri(D)], D[upper.tri(D)])
#> [1] 0.9757443
```

## Bayesian unfolding

[`bunfld()`](https://rootcoder007.github.io/crimlab/reference/bunfld.md):
The lognormal model of
[`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md)
for a respondent-by-stimulus matrix: delta\_{ij} N(x_i - z_j , 1 / tau),
N(0, 10^2) priors on every coordinate and tau U(0, 10) (Bakker and Poole
2013). Feeling thermometers T on 0-100 enter as delta = (100 - T) / 50.
Given the stimuli the respondents are independent, and the other way
round, so each block is slice sampled at once, dimension by dimension;
draws are aligned on the stimuli. The chain starts from the singular
value decomposition of the centred matrix. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

``` r

set.seed(1)
S <- matrix(runif(8, -0.5, 0.5), 4)
R <- matrix(runif(60, -0.5, 0.5), 30)
D <- sqrt(outer(R[, 1], S[, 1], "-")^2 + outer(R[, 2], S[, 2], "-")^2)
f <- bunfld(D, nsamp = 50, burn = 50)
cor(as.numeric(f$dist), as.numeric(D))
#> [1] 0.9766247
```

## Alpha-NOMINATE

[`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md):
Legislator i votes yea on roll call j with probability Phi(u\_{ij}),
u\_{ij} = Q\_{ij} + alpha (G\_{ij} - Q\_{ij}), where Q is the quadratic
and G the Gaussian utility difference of the distances to the yea and
nay outcomes, with weight 0.5 (Carroll et al. 2013): alpha = 1 is
Gaussian utility, alpha = 0 quadratic. beta has a flat prior on (0,
Inf), alpha on \[0, 1\]; coordinates carry the hierarchical penalty
-x’Sx/2 with S redrawn each sweep from an inverse Wishart. Every
coordinate, beta and alpha are slice sampled (Neal 2003: width 8, at
most 3 steps out); starts are uniform on (-1, 1), beta = 10, alpha =
0.7. Each draw is then centred on its legislators and rotated onto the
posterior mean, the `polarity` legislator positive. Roll calls whose
minority side is below `lop` and legislators with fewer than `minvotes`
votes are dropped first. Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

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

## Ordinal IRT (ordinal factor model)

[`ordirt()`](https://rootcoder007.github.io/crimlab/reference/ordirt.md):
Each ordinal response is a thresholded latent normal, y^**{ij} =
lambda*{j0} + lambda_j’ phi_i + e\_{ij}, y\_{ij} = c when gamma\_{j,c-1}
\< y^*\_{ij} \<= gamma\_{j,c}, with gamma\_{j1} = 0, phi_i N(0, I),
loadings lambda_j N(0, I / L_0) and a flat prior on the free cutpoints
(Quinn 2004). Each sweep updates the free cutpoints of every item by the
Metropolis-Hastings step of Cowles (1996) with y^\* integrated out (step
sizes tuned during the burn-in only), then draws y^\*, the loadings and
the ideal points from their full conditionals. The chain starts from the
leading principal components of the mean-filled responses. Each draw is
reflected so the first item loads positively (rotated onto the first
draw with several dimensions). Draws come from `seed`, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it.

``` r

set.seed(1)
th <- rnorm(60)
Y <- matrix(cut(outer(th, runif(6, 0.8, 1.5)) + rnorm(360),
                c(-Inf, -0.5, 0.5, Inf), labels = FALSE), 60, 6)
f <- ordirt(Y, nsamp = 50, burn = 50)
cor(f$ideal[, 1], th)
#> [1] 0.8995402
```

## Dynamic IRT (ideal points that drift)

[`dynirt()`](https://rootcoder007.github.io/crimlab/reference/dynirt.md):
The one-dimensional dynamic item response model of Martin and Quinn
(2002), made for courts: a vote is cast when z\_{jk} = -alpha_k + beta_k
theta\_{j,t(k)} + *{jk} \> 0; the ideal points follow a random walk,
theta*{j,0} N(e_0, E_0), theta\_{j,t} N(theta\_{j,t-1}, tau^2_j);
alpha_k N(a_0, 1/A_0), beta_k N(b_0, 1/B_0). Each sweep draws the latent
utilities, every case’s (alpha, beta), every justice’s path by forward
filtering and backward sampling and, when `c0` and `d0` are positive,
the evolution variances from IG((c_0 + T)/2, (d_0 + SS)/2). A draw is
reflected whenever the `anchor` member’s mean ideal point is negative.
The chain starts from the first principal component of the votes. Draws
come from `seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html)
takes it.

``` r

# 12 simulated justices over three terms of 15 cases
set.seed(7)
period <- rep(1:3, each = 15)
walk <- function(s) s + cumsum(c(0, rnorm(2, 0, 0.3)))
truth <- t(vapply(rnorm(12), walk, numeric(3)))
p <- pnorm(-rep(rnorm(45, 0, 0.5), each = 12) +
  rep(rnorm(45, 1.5, 0.5), each = 12) * truth[, period])
V <- matrix(rbinom(length(p), 1, p), 12)
f <- dynirt(V, period, nsamp = 100, burn = 50)
cor(as.numeric(f$theta), as.numeric(truth))
#> [1] 0.9013526
```

## Anchoring vignettes

[`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md):
Corrects self-ratings for differential item functioning with anchoring
vignettes, nonparametrically (King et al. 2004): each respondent’s
ratings of the vignettes, sorted, give that respondent’s thresholds
between categories, and the self-rating is recoded as the number of
thresholds at or below it. The spread of each threshold across
respondents measures the differential item functioning.

``` r

set.seed(1)
anchv(sample(1:5, 30, TRUE), matrix(sample(1:5, 90, TRUE), 30))$dif
#> [1] 1.0618786 1.0618786 1.2429109 0.8193072
```

## Bootstrap of scaling positions

[`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md):
Standard errors and percentile intervals of Aldrich-McKelvey stimulus
positions
([`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md))
or Basic Space weights
([`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md),
first dimension, of the matrix or of its transpose) by resampling
respondents with replacement (Efron and Tibshirani 1993). Draws come
from `seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html) takes
it. The sign of a Basic Space dimension is arbitrary, so its replicates
are best compared after aligning them.

``` r

set.seed(1)
Z <- matrix(rnorm(100), 20, 5) + rep(1:5, each = 20)
scboot(Z, nboot = 10)$se
#> [1] 0.09685089 0.13563335 0.07232775 0.09721725 0.07797385
```

## Normal vector of an outside measure

[`nrmvec()`](https://rootcoder007.github.io/crimlab/reference/nrmvec.md):
Projects an outside measure (a judge’s party, a legislator’s vote share)
into an ideal-point space: the least-squares regression of the measure
on the coordinates, its slopes scaled to unit length give the direction
along which the measure rises, with its angle and the R^2 (Armstrong et
al. 2021, section 2.6).

``` r

set.seed(1)
nrmvec(matrix(rnorm(20), 10), rnorm(10))$degrees
#> [1] 153.4453
```

## Cutting lines of votes in two dimensions

[`cutlns()`](https://rootcoder007.github.io/crimlab/reference/cutlns.md):
The ends of each vote’s cutting line n_1 x + n_2 y = c across the
horizontal range `xlim` (or from -10 to 10 vertically when the line is
vertical), its midpoint and the angle of its normal: the lines of a
Coombs mesh (Poole 2005). The plot itself is not drawn here.

``` r

set.seed(1)
cutlns(matrix(rnorm(6), 3), c(0.1, -0.2, 0))$ends
#>      [,1]        [,2] [,3]       [,4]
#> [1,]   -1 -0.33000699    1  0.4553768
#> [2,]   -1 -0.04963973    1 -1.1642922
#> [3,]   -1  1.01847753    1 -1.0184775
```

## Wordfish: positions from word counts

[`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md):
The Poisson scaling model of Slapin and Proksch (2008), y\_{ij}
((psi_i + alpha_j + beta_j omega_i)), for documents (judgments, briefs,
platforms) by their word counts: alternating Newton steps on the
document positions omega (standardised after each sweep, the word
parameters adjusted to keep the linear predictor) and on the word
parameters, from normal starts (seed 42, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) takes it), until the
positions change by a relative `tol`. The sign of the scale is
arbitrary.

``` r

set.seed(1)
wrdfsh(matrix(rpois(600, 5), 20), maxit = 20)$pos
#>  [1]  1.15570659  0.25111445  0.01650556 -1.49106129  0.13088144 -1.07365520  1.05617759
#>  [8] -0.83530154 -0.99936159  0.25353260 -0.31732379 -1.84826884 -0.20769336  1.06437931
#> [15]  1.98755859  0.68871590  0.01462927 -0.99868738  1.10180015  0.05035155
```
