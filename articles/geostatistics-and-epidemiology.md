# Geostatistics and spatial epidemiology

Variograms and kriging, and area-level models of rates with spatial
dependence.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Geostatistics

### Empirical semivariogram

[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md):
The classical estimator of Matheron (1963): half the mean squared
difference of the values over the pairs of points whose distance falls
in each of `nbins` equal bins of (0, cutoff\] (the first bin closed at
0). Bins without pairs are left out. The semivariogram describes how
dissimilar values become with distance, the first step of kriging: rates
of a crime, a disease or a survey outcome measured at places.

``` r

# simulated values with a trend from west to east
set.seed(1)
p <- cbind(x = runif(60, 0, 10), y = runif(60, 0, 10))
v <- p[, 1] * 0.5 + rnorm(60, 0, 0.3)
vgemp(p, v = v)
#>         dist      gamma np
#> 1  0.3883710 0.03235575 14
#> 2  0.6472850 0.10202814 16
#> 3  0.9061991 0.14740942 25
#> 4  1.1651131 0.17579384 27
#> 5  1.4240271 0.12605216 28
#> 6  1.6829411 0.20833025 43
#> 7  1.9418551 0.27639636 43
#> 8  2.2007691 0.37523216 57
#> 9  2.4596831 0.46234587 48
#> 10 2.7185972 0.54284954 60
#> 11 2.9775112 0.68058945 69
#> 12 3.2364252 0.78158842 69
#> 13 3.4953392 0.74873781 60
#> 14 3.7542532 0.84625342 76
```

### Fit a variogram model by Gaussian maximum likelihood

[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md):
Fits the nugget, partial sill and range of a spherical, exponential or
Gaussian covariance model by maximising the Gaussian log-likelihood of
the values about a constant mean (Mardia and Marshall 1984). The search
is the Nelder-Mead simplex of R’s
[`optim()`](https://rdrr.io/r/stats/optim.html) (at most 800
evaluations) on the logarithms of the three parameters, from 0.1 and 0.9
of the variance of the values and a quarter of the largest distance. The
covariance matrix gets `1e-8` of the variance on its diagonal so that it
stays positive definite.

``` r

set.seed(1)
p <- cbind(x = runif(40, 0, 10), y = runif(40, 0, 10))
v <- p[, 1] * 0.5 + rnorm(40, 0, 0.3)
vgfit(p, v = v)
#> vgfit: exponential model by maximum likelihood
#>       nugget        psill        range       loglik 
#>  0.006099213  8.719291632 71.862753310 11.408010219
```

### Ordinary kriging

[`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md):
Predicts the values at new locations `at` by ordinary kriging: the best
linear unbiased predictor under a variogram model with an unknown
constant mean (Cressie 1993, sec. 3.2), with its kriging variance. The
variogram is `vgm`, a list such as
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)
returns; without it, an exponential model is fitted by weighted least
squares to the 15-bin
[`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md)
semivariogram with the weights `np / h^2` of Cressie (1985), which is
fast but rougher than
[`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md).

``` r

set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- sin(3 * p[, 1]) + rnorm(60, sd = 0.1)
okrge(p, v = v, at = cbind(0.5, 0.5))
#>      pred        var
#> 1 1.06704 0.01030984
okrge(p, v = v, at = cbind(0.5, 0.5),
      vgm = list(model = "exponential", nugget = 0.01, psill = 0.3,
                 range = 0.4))
#>       pred        var
#> 1 1.062215 0.04918238
```

### Ordinary block kriging

[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md):
Predicts the mean value over each block (an area such as a district) by
ordinary kriging with point-to-block covariances averaged over the
block, C(s_i, B) = (1 / \|B\|) int_B C(s_i, u) du, under the exponential
covariance (sill - nugget) exp(-h / range), plus `nugget` at distance 0
(Schabenberger and Gotway 2005, ch. 4). A block is a box, given as a 2 x
2 matrix of its lower and upper corners (rows) and integrated on a
regular grid of `round(sqrt(nq))^2` points, or a matrix of the points
that discretise it.

``` r

set.seed(1)
p <- cbind(x = runif(30), y = runif(30))
v <- p[, 1] + rnorm(30, sd = 0.1)
bkrig(p, v = v, blocks = list(rbind(c(0, 0), c(0.5, 0.5)),
                              rbind(c(0.5, 0.5), c(1, 1))),
      range = 0.3)
#>         est        se
#> 1 0.3357770 0.1852905
#> 2 0.7159535 0.2114923
```

### Leave-one-out cross-validation of ordinary kriging

[`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md):
Predicts each value from all the others by ordinary kriging with the
exponential covariance of
[`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md),
and summarises the residuals by the mean squared prediction error MSPE =
(1 / n) sum (Z(s_i) - Z\_-i(s_i))^2, its root and the mean absolute
error (Schabenberger and Gotway 2005, ch. 4). Comparing it across
covariance parameters chooses among them.

``` r

set.seed(1)
p <- cbind(x = runif(30), y = runif(30))
krigcv(p, v = rnorm(30), range = 0.2)
#> krigcv: leave-one-out kriging of 30 values
#>      mspe     rmspe       mae 
#> 0.7067174 0.8406648 0.7098562
```

### Polynomial trend surface

[`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md):
Fits the mean of the values as a polynomial in the coordinates by
ordinary least squares, mu(s) = sum_k beta_k f_k(s), the `f_k` the
monomials up to degree `order`: `1`; `x`, `y`; `x2`, `y2`, `xy`; `x3`,
`y3`, `x2y`, `xy2` (Schabenberger and Gotway 2005, ch. 2). The large-
scale trend is often removed this way before a variogram is fitted to
the residuals.

``` r

set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
v <- 1 + 2 * p[, 1] - p[, 2]^2 + rnorm(50, sd = 0.1)
f <- trsurf(p, v = v)
f
#> trsurf: trend surface of order 2 on 50 points, R-squared 0.9818
#>           1           x           y          x2          y2          xy 
#>  1.09625518  2.11305475 -0.73267272 -0.08022062 -0.33750250  0.14443429
summary(f)
#>   term         est         se          z
#> 1    1  1.09625518 0.07968901 13.7566670
#> 2    x  2.11305475 0.21086371 10.0209505
#> 3    y -0.73267272 0.25326074 -2.8929581
#> 4   x2 -0.08022062 0.18258998 -0.4393484
#> 5   y2 -0.33750250 0.21097762 -1.5997076
#> 6   xy  0.14443429 0.18005051  0.8021877
```

### Spatial linear model by profile maximum likelihood

[`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md):
The Gaussian spatial linear model v = Z beta + delta + e, `delta` a
Gaussian process with covariance sigma2 exp(-h / phi), without a nugget
(Schabenberger and Gotway 2005, ch. 5). The range `phi` is found by
maximising the profile likelihood with Brent’s method on `log(phi)` over
\[log(max(h / 100, 0.001)), log(max(3 h, 1))\], `h` the largest
distance, as R’s [`optimize()`](https://rdrr.io/r/stats/optimize.html)
runs it; `beta` follows by generalised least squares. The standard
errors use `sigma2` estimated with `n - p` degrees of freedom.

``` r

set.seed(1)
p <- cbind(x = runif(40), y = runif(40))
z <- cbind(1, p[, 1])
sglmx(1 + 2 * p[, 1] + rnorm(40, sd = 0.2), z, p)
#> sglmx: spatial linear model, exponential covariance, 40 points
#>         z1         z2     sigma2        phi 
#> 0.97349153 2.09618531 0.03478370 0.01495528
```

### Spatial linear mixed model by REML

[`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md):
The model v = Z beta + delta + e, delta ~ N(0, sigma2 R_phi) with the
exponential correlation `R_phi = exp(-h / phi)` and a nugget e ~ N(0,
tau2 I), tau2 = nu sigma2, by restricted maximum likelihood (Patterson
and Thompson 1971; Schabenberger and Gotway 2005, ch. 5). The search is
R’s Nelder-Mead ([`optim()`](https://rdrr.io/r/stats/optim.html), 400
evaluations, relative tolerance `1e-6`) on (log phi, log nu) from (log(h
/ 3), log(0.1)), `h` the largest distance (at least 1).

``` r

set.seed(1)
p <- cbind(x = runif(40), y = runif(40))
smixm(1 + 2 * p[, 1] + rnorm(40, sd = 0.3), cbind(1, p[, 1]), p)
#> smixm: spatial linear mixed model by REML, 40 points
#>         z1         z2     sigma2       tau2        phi 
#> 0.96148075 2.14380681 0.06465755 0.01544588 0.03247176
```

### A spatial data object

[`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md):
Holds points with their attributes and their coordinate reference
system: the coordinates come from two columns of a data frame, the other
columns are kept as attributes. The system is an EPSG code, never a PROJ
string: 4326 (WGS 84 longitude and latitude, in degrees), 3857 (Web
Mercator, in metres) or another code, taken to be in metres.
[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md)
returns the coordinates,
[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md)
the system,
[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md)
the specification of a web map and
[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md)
a sample of the points.

``` r

d <- data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7), v = c(1, 2))
s <- spobj(d, c("lon", "lat"))
s
#> spobj: 2 points, EPSG:4326 (degrees)
format(s)
#> [1] "spobj: 2 points, EPSG:4326 (degrees)"
summary(s)
#>   coord   min   max
#> 1   lon -79.4 -79.3
#> 2   lat  43.6  43.7
```

### The coordinates of a spatial data object

[`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md):
The coordinates of a spatial data object

``` r

spxy(spobj(data.frame(lon = 1:3, lat = 1:3)))
#>      lon lat
#> [1,]   1   1
#> [2,]   2   2
#> [3,]   3   3
```

### The coordinate reference system of a spatial data object

[`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md):
Gives the system as its EPSG code and a label in the style of WKT2,
never as a PROJ string.

``` r

crsof(spobj(data.frame(lon = 1, lat = 1)))
#> $epsg
#> [1] 4326
#> 
#> $wkt
#> [1] "GEOGCRS[\"WGS 84\", EPSG:4326]"
```

### The specification of a web map

[`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md):
Describes an interactive map of the points as data, for a mapping
library to draw: the tile provider, the centre (the mean longitude and
latitude), a zoom level and the points in longitude and latitude. Web
Mercator coordinates (EPSG:3857) are taken back to longitude and
latitude by the inverse of the spherical Mercator projection.

``` r

lfspec(spobj(data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7))))
#> $provider
#> [1] "OpenStreetMap"
#> 
#> $interactive
#> [1] TRUE
#> 
#> $center
#>    lon    lat 
#> -79.35  43.65 
#> 
#> $zoom
#> [1] 8
#> 
#> $points
#>        lon  lat
#> [1,] -79.4 43.6
#> [2,] -79.3 43.7
```

### Sample points, by local density or uniformly

[`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md):
Draws `size` of the points without replacement. With `dens = TRUE` each
draw picks a point with probability proportional to its local density
sum_j exp(-(d_ij / b)^2), `b` the median distance between points
(great-circle distance for EPSG:4326), so that dense areas are sampled
more; otherwise every point is equally likely.

``` r

s <- spobj(expand.grid(lon = 1:5, lat = 1:5))
spsamp(s, 5, seed = 1)
#> [1]  8 20 10  6 18
#> attr(,"prob")
#>  [1] 0.02813274 0.03613227 0.03913468 0.03613227 0.02813274 0.03614678 0.04642113
#>  [8] 0.05027751 0.04642113 0.03614678 0.03916871 0.05029579 0.05447252 0.05029579
#> [15] 0.03916871 0.03618351 0.04645522 0.05031121 0.04645522 0.03618351 0.02818869
#> [22] 0.03618420 0.03918600 0.03618420 0.02818869
```

## Spatial epidemiology

### Bayes outbreak detection

[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md):
The Bayes detector of `surveillance::algo.bayes()` (Riebler 2004; Hoehle
2007) for a series of counts: at each time `t` the reference counts are
the `w` before it (`acty = TRUE`) and, in each of the `b` years before,
the 2 w + 1 counts around the same time of year (`freq` counts per
year); missing counts are left out. With `s` their sum and `n` their
number, the posterior predictive distribution of the count at `t` under
a Poisson model with a flat prior is negative binomial of size `s + 0.5`
and probability `n / (n + 1)`; its `1 - alpha` quantile is the upper
bound, and a count above it raises an alarm. The same series may be
calls for service, emergency visits or reported crimes of one kind.

``` r

bayout(c(3, 5, 2, 4, 6, 3, 4, 12), w = 6, at = 8)
#>   t ub alarm
#> 1 8  8  TRUE
```

### Exposure within a buffer

[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md):
For each receptor (a home, a school), the number of sources (bars,
polluting sites, hot spots) within distance `r`, their summed weights
and the distance to the nearest: the buffer measure of exposure
assessment. Distances are planar, or great-circle metres when
`lonlat = TRUE` (the coordinates then longitude and latitude in
degrees).

``` r

rec <- cbind(x = c(0, 1, 5), y = c(0, 1, 5))
src <- cbind(x = c(0.5, 4), y = c(0.5, 4))
bufexp(rec, src = src, r = 1.5)
#>   count wsum      near
#> 1     1    1 0.7071068
#> 2     1    1 0.7071068
#> 3     1    1 1.4142136
```

### Exposure by kernel

[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md):
The kernel-weighted density of sources at each receptor, sum_j w_j
K(d_ij / h) / h^2, with the bivariate normal kernel exp(-u^2 / 2) / (2
pi) or the quartic kernel `3 / pi (1 - u^2)^2` (for `u < 1`, 0 beyond)
of Bailey and Gatrell (1995). Planar coordinates.

``` r

kerexp(rbind(c(0.2, 0.3)), src = rbind(c(0.25, 0.35), c(0.6, 0.8)),
       h = 0.5, kern = "quartic")
#> [1] 3.668458
```

### The Leroux precision matrix

[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md):
The precision of the area effects in the model of Leroux, Lei and
Breslow (2000), tau (rho Q + (1 - rho) I), with `Q = D - W` the
Laplacian of the neighbour weights: independent effects at `rho = 0`,
the intrinsic conditional autoregression at `rho = 1`.

``` r

lerprc(rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0)), 0.5, 2)
#>      [,1] [,2] [,3]
#> [1,]    2   -1    0
#> [2,]   -1    3   -1
#> [3,]    0   -1    2
```

### The BYM2 area structure

[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md):
The scaled structure of the BYM2 model of Riebler et al. (2016): the
generalised inverse `G` of the Laplacian `Q` of the neighbour weights
(from its eigenvalues above `tol` times the largest), the scaling factor
`s`, the geometric mean of the diagonal of `G`, so that the structured
effect has a generalised variance of 1, and the covariance of the area
effects ((1 - phi) I + phi G / s) / tau. The graph must be connected and
every area must have a neighbour.

``` r

w <- rbind(c(0, 1, 1, 0), c(1, 0, 1, 0), c(1, 1, 0, 1), c(0, 0, 1, 0))
bym2sc(w, 2, 0.4)$scale
#> [1] 0.3565926
```

### Ecological regression with expected counts

[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md):
Regression of area counts `y` on covariates with the offset `log(e)`,
`e` the expected counts (from indirect standardisation, say), so that
the coefficients act on the standardised ratio: Poisson by iteratively
reweighted least squares; negative binomial, alternating Newton steps on
the score of `theta` with reweighted fits of the coefficients; or
zero-inflated Poisson (Lambert 1992) by EM, the zero part logistic on an
intercept and the columns of `z`.

``` r

# simulated counts and expected counts
set.seed(1)
x <- rnorm(30)
e <- runif(30, 3, 8)
y <- rpois(30, e * exp(0.2 + 0.4 * x))
ecolrg(y, e, x)
#> ecolrg: 30 counts, log-likelihood -67.3106
#> [1] 0.1404659 0.4851843
```

### Eigenvector spatial filtering

[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md):
Adds Moran eigenvectors of the neighbour weights to the regression of
`v` on an intercept and `x` by forward selection, so that the residuals
lose their spatial autocorrelation (Griffith 2003; Tiefelsdorf and
Griffith 2007). The Moran eigenvectors are the eigenvectors of the
symmetrised weights (W + W’) / 2 orthogonal to the constant (computed in
the Helmert basis), each with Moran’s I n / S0 lambda; the candidates
are those with Moran’s I above `-1 / (n - 1)` (`"positive"`, the usual
choice), below it, or all. Selection by `method`: `"moran"` takes the
vector that most reduces the absolute residual Moran’s I until it is
below `tol`; `"aic"` and `"press"` (the leave-one-out prediction error)
while they fall; `"r2"` while R-squared rises by `tol` or more. At least
two residual degrees of freedom are kept.

``` r

w <- 1 * (abs(outer(1:8, 1:8, "-")) == 1)
v <- c(1, 2, 2.5, 4, 3.5, 3, 1.5, 1)
evfilt(v, w, method = "aic")
#> evfilt: 1 eigenvectors selected, residual Moran's I -0.7032, R-squared 0.9199
```

### The metapopulation SIR model

[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md):
The deterministic SIR model on patches (Keeling and Rohani 2008,
sec. 7.3): the force of infection on patch `i` is beta sum_j C_ij I_j /
N_j with the coupling C = (1 - kappa) I + kappa W~, W~ the
row-standardised neighbour weights, so that a share `kappa` of the
pressure comes from neighbours; recovery is at rate `gamma`. Integrated
by the classical fourth-order Runge-Kutta method in steps of `dt`. The
same equations describe the spread of an outbreak of violence between
neighbourhoods in contagion models.

``` r

w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
spsir(c(10, 0, 0), w, rep(1000, 3), beta = 0.6, gamma = 0.2,
      tmax = 60, dt = 0.5)
#> spsir: 3 patches to time 60, attack rate 0.9399, R0 3
```

### Spatial frailty and cure models

[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`sfcure()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md):
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md)
fits the Weibull proportional hazards model with log-frailties of the
areas under an intrinsic conditional autoregression, h(t) = rho
t^(rho - 1) exp(x beta + u_area) (Banerjee, Wall and Carlin 2003);
[`sfcure()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md)
the promotion-time cure model (Chen, Ibrahim and Sinha 1999),
`S(t) = exp(-theta (1 - exp(-G(t))))`, theta = exp(x beta + u_area) and
G(t) = exp(c + rho log t), whose cure fraction is `exp(-theta)`. Both
maximise the penalised likelihood, the penalty tau u’Qu / 2 with `Q` the
Laplacian of the neighbour weights, by Newton-Raphson from 0 on the
analytic gradient with a Hessian by central differences, the step halved
while the objective falls, until a step is below `tol` (200 iterations
at most). The standard errors come from the last Hessian; the first is
that of the log shape.

``` r

# simulated survival times in three areas
w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
i <- 0:29
x <- sin(i)
t <- 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * x - 0.3 * (i %% 3))
f <- sfrail(t, as.numeric(i %% 5 != 0), x, i %% 3 + 1, w)
f
#> sfrail: Weibull frailty model, shape 3.065, penalised log-likelihood -20.2083
#> [1]  0.9000010 -1.9361730 -1.2253279 -0.8237804
summary(f)
#>       term        est        se         z
#> 1 logshape  1.1201896 0.1648807  6.793941
#> 2       x1  0.9000010 0.3560158  2.527981
#> 3   frail1 -1.9361730 0.4750970 -4.075321
#> 4   frail2 -1.2253279 0.4050962 -3.024783
#> 5   frail3 -0.8237804 0.3803217 -2.166009
```

### The spatial encoding of graph attention

[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md):
Turns a matrix of hop distances between the nodes of a network (areas,
places, people) into the bias that Graphormer (Ying et al. 2021) adds to
the attention logits: a learnable value per distance, `b[d + 1]` (the
last for every larger distance), and `ub` for pairs that cannot reach
each other (distance `-1`). A bias inside the softmax keeps distant
nodes reachable but discouraged.

``` r

d <- rbind(c(0, 1, 2), c(1, 0, -1), c(2, -1, 0))
spbias(d, b = c(1, 0.5, 0.25))
#> $bias
#>      [,1]  [,2]   [,3]
#> [1,] 1.00   0.5   0.25
#> [2,] 0.50   1.0 -10.00
#> [3,] 0.25 -10.0   1.00
#> 
#> $ub
#> [1] -10
```
