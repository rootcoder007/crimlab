# Spatial geometry and centrography

Distances, projections, centres, ellipses and the other geometry that
spatial analysis of crime starts from.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Great-circle distance

[`hvdist()`](https://rootcoder007.github.io/crimlab/reference/hvdist.md):
The distance along the sphere between pairs of points given in longitude
and latitude. The central angle is computed as the `atan2` of the cross
and dot products of the two unit vectors, written so that nothing
cancels: accurate to a few units in the last place from coincident
points to antipodes (the spherical case of Vincenty 1975). The Earth is
taken as a sphere of mean radius `r`; distances on the ellipsoid differ
by up to about 0.5 per cent.

``` r

# Toronto to Montreal, metres
hvdist(c(-79.38, 43.65), c(-73.57, 45.50))
#> [1] 503922
pts <- data.frame(lon = c(0, 90, 180), lat = c(0, 0, 0))
hvdist(pts, c(0, 0))
#> [1]        0 10007557 20015114
```

## Equal-area projection of longitude and latitude, and its inverse

[`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md),
[`eqinv()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md):
[`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
maps longitude and latitude to planar metres by the Lambert azimuthal
equal-area projection of the sphere about a centre;
[`eqinv()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
maps them back. Areas are exact, so intensities, nearest-neighbour
expectations and K functions computed on the result are right; lengths
are scaled by `cos(c / 2)` along the radius from the centre and its
inverse across it (`c` the central angle), so within 100 km of the
centre every length is right to 3.1e-5 relative. Every planar function
of the package takes the projected coordinates; those with a `lonlat`
argument call this projection themselves.

``` r

p <- eqproj(c(-79.40, -79.38, -79.35), c(43.64, 43.66, 43.65))
p
#>            x             y
#> 1 -1877.6518 -1111.6869087
#> 2  -268.1467  1111.9561852
#> 3  2145.5306     0.3446339
a <- attr(p, "centre")
eqinv(p, lon0 = a[1], lat0 = a[2])
#>      lon   lat
#> 1 -79.40 43.64
#> 2 -79.38 43.66
#> 3 -79.35 43.65
```

## Mean centre and dispersion of points

[`cgstat()`](https://rootcoder007.github.io/crimlab/reference/cgstat.md):
The (weighted) mean centre of a point pattern, the standard deviations
of the two coordinates, the standard distance deviation (the root mean
squared distance from the mean centre), and the geometric and harmonic
means of the coordinates (Bachi 1959; Neft 1966). The mean centre is the
point that minimises the sum of squared distances to the points. The
dispersions divide by the total weight less `k`: `k = 1` makes the
squared values unbiased for the variances; `k = 2` is the convention of
some crime-mapping texts; `k = 0` divides by the total weight. That the
mean centre minimises the weighted sum of squared distances is proved in
Lean: see `proof("cgstat")`.

``` r

set.seed(1)
p <- cbind(x = rnorm(50, 100, 10), y = rnorm(50, 200, 5))
cgstat(p)
#>          w         mx         my         sx         sy        sdd         gx         gy 
#>  50.000000 101.004483 200.586632   8.313939   4.844140   9.622228 100.653418 200.529543 
#>         hx         hy 
#> 100.284986 200.472685
cgstat(p[, 1], p[, 2], w = rep(1:2, 25))
#>          w         mx         my         sx         sy        sdd         gx         gy 
#>  75.000000 100.685776 200.457192   8.553285   4.753294   9.785321 100.307901 200.401917 
#>         hx         hy 
#>  99.908680 200.346969
cgstat(p, k = 0)  # dispersions over n, not n - 1
#>          w         mx         my         sx         sy        sdd         gx         gy 
#>  50.000000 101.004483 200.586632   8.230380   4.795453   9.525520 100.653418 200.529543 
#>         hx         hy 
#> 100.284986 200.472685
```

## Median centre of points

[`cgmedn()`](https://rootcoder007.github.io/crimlab/reference/cgmedn.md):
The median of the `x` coordinates and the median of the `y` coordinates
(the mean of the two middle values for an even number of points). It
depends on the orientation of the axes: rotating the points does not
rotate it.
[`webpt()`](https://rootcoder007.github.io/crimlab/reference/webpt.md)
gives the rotation-invariant centre of minimum distance.

``` r

cgmedn(c(1, 2, 3, 10), c(5, 1, 2, 2))
#>   x   y 
#> 2.5 2.0
```

## Centre of minimum distance

[`webpt()`](https://rootcoder007.github.io/crimlab/reference/webpt.md):
The point `c` that minimises the weighted sum of distances sum(w \*
\|p - c\|) to the points: the Weber point or spatial median, which,
unlike the median centre, does not depend on the axes. It is found by
Weiszfeld’s iteration with the modification of Vardi and Zhang (2000),
which keeps every step defined and the sum non-increasing when an
iterate lands on a data point. The optimality gap certifies the answer:
it is 0 exactly at a minimiser. The minimiser is unique unless all the
points lie on one line with the weight split evenly.

``` r

webpt(c(0, 1, 0, 10), c(0, 0, 1, 10))
#>            x            y          obj          gap         iter         code 
#> 5.000000e-01 5.000000e-01 1.555635e+01 4.086973e-10 2.200000e+01 1.000000e+00
```

## Standard deviational ellipse

[`sdell()`](https://rootcoder007.github.io/crimlab/reference/sdell.md):
The ellipse whose axes lie along the eigenvectors of the weighted
covariance matrix of the points about their mean centre, with semi-axes
`c` times the square roots of its eigenvalues (Wang, Shi and Miao 2015).
It is the ellipse of the rotated-axes formulas of Lefever (1926) and
Ebdon (1988) once the same divisor is used; texts differ only in the
divisor (`k`) and in the multiplier (some multiply by `sqrt(2)`). The
curve of the standard deviation by direction that Lefever drew is not
this ellipse (Gong 2002); the ellipse is returned. For bivariate normal
points the ellipse with `c = sqrt(-2 * log(1 - p))` covers a share `p`
of them: `c = 1` covers 39.3 per cent, `c = 2` 86.5 per cent. The
parallel-axis identity behind the ellipse’s centre is proved in Lean:
see `proof("sdell")`.

``` r

set.seed(2)
x <- rnorm(200)
sdell(x, 0.5 * x + rnorm(200, sd = 0.3), c = 2)
#>            mx            my             a             b         angle          area 
#> -0.0007417526  0.0681841044  2.4146550213  0.5262845280  0.4870207069  3.9923220528
```

## Convex hull of points

[`cvhull()`](https://rootcoder007.github.io/crimlab/reference/cvhull.md):
The smallest convex polygon holding every point, by Andrew’s monotone
chain with an exact orientation predicate (Shewchuk 1997), so
near-collinear points never make a wrong turn. Points on an edge are not
vertices; of coincident points the first is reported.

``` r

cvhull(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5))
#> $idx
#> [1] 1 2 3 4
#> 
#> $area
#> [1] 1
#> 
#> $perim
#> [1] 4
```

## Nearest-neighbour index

[`nnidx()`](https://rootcoder007.github.io/crimlab/reference/nnidx.md):
The mean distance from each point to its k-th nearest neighbour, against
its expectation under complete spatial randomness (CSR) in the rectangle
`win`: the nearest-neighbour index of Clark and Evans (1954) and its
k-order extension (Thompson 1956), with Donnelly’s (1978) edge-corrected
expectation for the first order and, when `nsim > 0`, a Monte Carlo test
against `nsim` patterns of as many uniform points in the rectangle. The
Monte Carlo p-values have exact size; the z statistics ignore that
nearest-neighbour distances are dependent and, uncorrected, the edge.
The index tests only the k-th neighbour distances, not complete spatial
randomness as a whole.

``` r

set.seed(1)
x <- runif(200)
y <- runif(200)
nnidx(x, y, win = c(0, 1, 0, 1), kmax = 3, nsim = 99)
#> $order
#>   k       mean        exp    ratio pclu preg         sd
#> 1 1 0.03753432 0.03535534 1.061631 0.84 0.17 0.01848106
#> 2 2 0.05462585 0.05303301 1.030035 0.44 0.57 0.01925094
#> 3 3 0.06778844 0.06629126 1.022585 0.19 0.82 0.01949659
#> 
#> $ce
#>  donnelly     ratio         z        zc 
#> 0.0364416 1.0299854 1.6674068 0.8361710 
#> 
#> $dist
#> NULL
```

## Ripley’s K function

[`ripk()`](https://rootcoder007.github.io/crimlab/reference/ripk.md):
Ripley’s K function of a point pattern in the rectangle `win` at radii
`r`, and Besag’s `L = sqrt(K / pi)`. Under complete spatial randomness
K(r) = pi r^2 and `L(r) = r`. Four estimators: `"none"` (no edge
correction), `"border"` (only points at least `r` from the edge are
centres), `"modif"` (the modified border estimator, with the area of the
eroded window; `NaN` once `r` reaches half a side) and `"iso"` (Ripley’s
isotropic correction, each pair weighted by the inverse of the share of
the circle through the other point that lies in the rectangle, computed
exactly). The normalisation is `n (n - 1)`, as in spatstat.

``` r

set.seed(1)
ripk(runif(100), runif(100), win = c(0, 1, 0, 1), r = seq(0, 0.25, 0.05))
#>      r           k          l
#> 1 0.00 0.000000000 0.00000000
#> 2 0.05 0.007686885 0.04946525
#> 3 0.10 0.033207387 0.10281167
#> 4 0.15 0.080766297 0.16033936
#> 5 0.20 0.131964756 0.20495289
#> 6 0.25 0.201545381 0.25328618
```

## Envelopes of Ripley’s L function and a test of complete spatial randomness

[`ripenv()`](https://rootcoder007.github.io/crimlab/reference/ripenv.md):
`L(r) - r` for the data, its pointwise minimum, maximum and mean over
`nsim` patterns of as many uniform points in the rectangle, and a global
Monte Carlo test on the maximum absolute deviation T = max_r \|L(r) -
r\|, whose p-value (1 + \#{T_sim \>= T}) / (nsim + 1) has exact size
(Ripley 1977; Diggle 2013). The pointwise envelopes are not a test of
the whole curve: each radius alone falls outside them with probability
`2 / (nsim + 1)` under randomness.

``` r

set.seed(1)
ripenv(runif(100), runif(100), win = c(0, 1, 0, 1),
       r = seq(0.01, 0.25, 0.02), nsim = 39)$test
#>          t          p 
#> 0.01033936 0.15000000
```

## Fuzzy mode

[`fzmode()`](https://rootcoder007.github.io/crimlab/reference/fzmode.md):
For each point, the number of points (itself and coincident ones
included) within distance `r` of it. With `r = 0` it counts coincident
points, the mode of exact locations; with `r > 0` it is the fuzzy mode
of the crime-mapping literature, a fixed-radius scan centred on each
incident. Counts of nearby points overlap by design: each point counts
its own neighbourhood.

``` r

fzmode(c(0, 0, 1, 5), c(0, 0, 0, 5), r = 1)
#> [1] 3 3 3 1
```

## K-means partition of points

[`kmpart()`](https://rootcoder007.github.io/crimlab/reference/kmpart.md):
Partitions the points into `k` groups by Lloyd’s iterations: each point
joins its nearest centre (the lowest-numbered on ties), then each centre
moves to its points’ mean (`upd = "mean"`, K-means, minimising the
within-group sum of squares) or to their centre of minimum distance
(`upd = "weber"`, a spatial K-medians minimising the sum of distances),
until no point changes group or after `maxit` moves. An emptied group
keeps its centre. The result depends on the start: k-means++ seeding
(Arthur and Vassilvitskii 2007), the centres of the busiest grid cells
at least a given distance apart, or centres given.

``` r

set.seed(1)
x <- c(rnorm(30, 0), rnorm(30, 5))
y <- c(rnorm(30, 0), rnorm(30, 5))
kmpart(x, y, k = 2, seed = 1)$cen
#>            x        y  n       ss       sd
#> 1 5.13277458 5.113333 30 40.91834 31.93770
#> 2 0.08245817 0.110278 30 51.51169 33.83318
```

## Nearest-neighbour hierarchical clustering

[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md):
Groups points that lie closer together than a threshold distance, then
groups the groups, order by order (after Johnson 1967, as used for crime
hot spots). At the first order two points are linked when their distance
is below the threshold: `fixd`, or the mean nearest-neighbour distance
expected under complete spatial randomness in an area `area` plus `z`
standard errors, 0.5 sqrt(A / m) + z sqrt((4 - pi) A / (4 pi)) / m for
`m` units, with `z = qnorm(p)` (`p = 0.5` gives the expectation itself,
smaller `p` a stricter threshold). Points with a linked neighbour are
taken in decreasing number of links; each one not yet clustered takes
every unclustered point linked to it; groups of at least `minpts` points
are kept. Their centres of minimum distance become the units of the next
order, with the threshold recomputed for their number, while at least
four remain, up to `maxord` orders. The iterative reassignment some
descriptions add after the first pass is not done: its rule is not
stated unambiguously.

``` r

set.seed(1)
x <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
y <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
nnhclu(x, y, area = 225, p = 0.05, minpts = 5)$clusters
#>   order points units         mx           my         wx          wy         a         b
#> 1     1     40    40 0.03845448 -0.003796574 0.06246601 -0.03941408 0.2470327 0.2266867
#> 2     1     36    36 6.00748914  5.922730801 5.99404607  5.93418792 0.3063012 0.2568727
#> 3     1      5     5 6.13610851  6.711097797 6.09427553  6.66090240 0.1958946 0.0298542
#>      angle    thresh
#> 1 2.820476 0.6309155
#> 2 2.060319 0.6309155
#> 3 0.769340 0.6309155
```

## Search-circle hot-spot areas

[`stacel()`](https://rootcoder007.github.io/crimlab/reference/stacel.md):
Hot-spot areas from search circles (Block and Block 1995): circles of
radius `r` are centred on the nodes of an `nx` by `ny` grid spanning the
rectangle `win` (corners included), the points in each are counted, the
circles holding at least two points are ranked by count and the first
`top` are kept; kept circles that share a point are merged, repeatedly,
and a hot-spot area is the set of points in a merged group, kept if it
has at least `minpts`. Descriptions differ on the circle’s radius and on
the grid; both are the caller’s choice here.

``` r

set.seed(1)
x <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
y <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
stacel(x, y, win = c(0, 10, 0, 10), r = 0.75, top = 5, minpts = 5)$clusters
#>   points       mx       my       wx       wy         a         b     angle
#> 1     35 2.027899 2.025648 2.055272 2.024903 0.3233517 0.2532551 0.1950933
```

## Risk-adjusted nearest-neighbour hierarchical clustering

[`rnnhcl()`](https://rootcoder007.github.io/crimlab/reference/rnnhcl.md):
Nearest-neighbour clustering (as
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md),
first order) with the threshold distance adjusted for a baseline
population at risk, so that incidents dense only because people are
dense there are not called a cluster. In
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)
the threshold is the mean nearest-neighbour distance under complete
spatial randomness plus `z` standard errors, 0.5 sqrt(A / n) + z
sqrt((4 - pi) A / (4 pi)) / n, which is `f / sqrt(lambda)` with
`lambda = n / A` and f = 0.5 + z sqrt((4 - pi) / (4 pi n)). Here the
uniform intensity is replaced by the intensity the `n` incidents would
have if they followed the baseline: lambda_i = n b(p_i) / B, where `b`
is the normal-kernel density (bandwidth `h`) of the baseline units with
sizes `bw` and `B` is their total, and two incidents `i` and `j` are
linked when \|p_i - p_j\| \< f / sqrt((lambda_i + lambda_j) / 2). Taking
the mean of the two intensities is this function’s reading of “the
threshold at the pair’s baseline density” (the published description
does not say how the two ends combine); a uniform baseline gives back
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)’s
threshold. Grouping is
[`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)’s
first order; higher orders are not formed.

``` r

set.seed(1)
town <- cbind(x = rnorm(400, 0, 1), y = rnorm(400, 0, 1))
inc <- rbind(town[1:60, ], cbind(x = rnorm(15, 3, 0.1),
                                 y = rnorm(15, 3, 0.1)))
r <- rnnhcl(inc, base = town, h = 0.5, minpts = 5)
r$clusters
#>   points       mx      my      wx       wy         a         b    angle
#> 1     15 2.973022 3.00705 2.97886 3.000247 0.1228748 0.1070088 1.886728
```

## Head-bang smoothing

[`hdbang()`](https://rootcoder007.github.io/crimlab/reference/hdbang.md):
Smooths values at points or areal units (rates by area centroid, say) by
Hansen’s (1991) head-bang rule, which removes isolated peaks and troughs
while keeping edges between regions of different level. For each unit
`i` the rule takes its `nn` nearest neighbours and every pair (j, k) of
them whose angle at `i` is at least `ang` degrees, so that `i` lies
roughly between `j` and `k` (a triple); with `ntrip > 0` only the
`ntrip` triples whose ends are nearest (smallest \|p_j - p_i\| + \|p_k -
p_i\|) are kept. Each triple gives a low value `min(v_j, v_k)` and a
high value `max(v_j, v_k)`; the low screen is the median of the lows,
the high screen the median of the highs, and the new value is the median
of the low screen, `v_i` and the high screen, so `v_i` moves only if it
lies outside both. With weights `w` (populations at risk, say) the
screens are weighted medians, each end weighted by its own unit’s
weight, after Mungiole, Pickle and Simonson (1999): the smallest value
whose cumulative weight reaches half the total. A unit with no triple
(on the outer edge of the points) keeps its value; Hansen’s extrapolated
end triples are not formed. Each of `iter` passes starts from the
previous pass’s values.

``` r

g <- expand.grid(x = 1:6, y = 1:6)
v <- rep(1, 36)
v[15] <- 50
hdbang(g, v = v)[15]
#> [1] 1
```

## Correlated random walk of a series of events

[`crwalk()`](https://rootcoder007.github.io/crimlab/reference/crwalk.md):
Treats the successive locations of a series (incidents attributed to one
offender, in time order) as a correlated random walk (Kareiva and
Shigesada 1983): steps of independent lengths and changes of heading
(turning angles) independent of each other and of the steps. The fit
gives the mean and mean squared step and the mean cosine and sine of the
turns; the expected squared displacement after `k` moves is Kareiva and
Shigesada’s for symmetric turning (mean sine 0), k E l^2 + 2 (E l)^2 c /
(1 - c) (k - (1 - c^k) / (1 - c)), next to the observed mean over every
start. The walk is then continued `nsim` times for `h` steps from the
last location and heading, each turn and step drawn independently and
uniformly from the observed ones (an empirical correlated random walk),
and the simulated positions summarise where the next events may fall. It
describes movement; it does not identify anyone, and a series with few
events gives a rough walk.

``` r

set.seed(1)
a <- cumsum(rnorm(20, 0, 0.5))
p <- cbind(x = cumsum(cos(a)), y = cumsum(sin(a)))
r <- crwalk(p, h = 3, nsim = 200, seed = 1)
r$fit
#>      steps         ml        ml2          c          s       head 
#> 19.0000000  1.0000000  1.0000000  0.8935287  0.1187582  1.9052388
r$msd
#>   k      obs      exp        x        y       sdx       sdy
#> 1 1 1.000000 1.000000 12.12662 10.09144 0.3882481 0.1794948
#> 2 2 3.787057 3.787057 11.60345 10.74197 0.7778394 0.4737584
#> 3 3 8.165410 8.170902 11.07326 11.25329 1.1794886 0.8584427
```

## Mantel’s space-time test

[`mantst()`](https://rootcoder007.github.io/crimlab/reference/mantst.md):
The Pearson correlation, over all pairs of events, between their spatial
dissimilarity (the distance, or `1 / (distance + cs)` when `cs > 0`,
Mantel’s reciprocal transform) and their temporal one (\|t_i - t_j\|, or
1 / (\|t_i - t_j\| + ct)), with a permutation test that permutes the
times over the fixed locations (Mantel 1967). Drawing distances and
intervals independently instead, as some descriptions do, does not keep
the dependence among the pairs of one distance matrix and gives no test
of stated size. Memory grows with the number of events, not pairs; time
with pairs times permutations.

``` r

set.seed(1)
x <- runif(60)
y <- runif(60)
mantst(x, y, t = x + rnorm(60, sd = 0.1), nsim = 199, seed = 1)
#>           r       pairs         pup         plo          p2 
#>    0.574889 1770.000000    0.005000    1.000000    0.005000
```

## Spatial-temporal moving average

[`stmavg()`](https://rootcoder007.github.io/crimlab/reference/stmavg.md):
For each event, the mean centre of the `span` events centred on it in
time order (ties in input order), the window shortened at the ends of
the series. It traces how the centre of a series of incidents moves over
time. It is a descriptive device of the crime-mapping literature with no
original paper to cite; the mean centre is that of
[`cgstat()`](https://rootcoder007.github.io/crimlab/reference/cgstat.md).

``` r

stmavg(x = 1:7, y = c(1, 1, 2, 3, 5, 8, 13), t = 7:1, span = 3)
#>    mx        my n
#> 1 1.5  1.000000 2
#> 2 2.0  1.333333 3
#> 3 3.0  2.000000 3
#> 4 4.0  3.333333 3
#> 5 5.0  5.333333 3
#> 6 6.0  8.666667 3
#> 7 6.5 10.500000 2
```

## Distance-decay functions for journey to crime

[`jtdecay()`](https://rootcoder007.github.io/crimlab/reference/jtdecay.md):
The decay of the chance of offending with distance from an offender’s
base, in the forms of Brantingham and Brantingham (1981), Rossmo (2000)
and Levine (2009): `"linear"`, f = max(0, A + B d), `par = c(A, B)`;
`"negexp"`, f = A exp(-b d), `par = c(A, b)` with `b >= 0`; `"normal"`,
f = A dnorm(d, mu, s), `par = c(A, mu, s)`; `"lognormal"`, f = A
dlnorm(d, mu, s), `par = c(A, mu, s)` (the lognormal density of `d`;
some texts use log d^2 and `1 / d^2`, which is not a density);
`"trunc"`, the truncated negative exponential, f = P d / dp up to the
peak `dp` and P exp(-C (d - dp)) beyond, `par = c(P, dp, C)`; and
`"table"`, values `tab` at distances 0, delta, 2 delta, …, linear in
between and 0 beyond, `par = delta` (as
[`jtcal()`](https://rootcoder007.github.io/crimlab/reference/jtcal.md)
tabulates).

``` r

jtdecay(c(0, 0.5, 1, 2, 5), type = "trunc", par = c(1, 0.5, 0.8))
#> [1] 0.00000000 1.00000000 0.67032005 0.30119421 0.02732372
```

## Journey-to-crime surface

[`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md):
Ranks the cells of a grid as candidate bases of an offender from the
locations of a series of incidents. `mode = "sum"` gives the classic
aggregate surface, the sum over incidents of the decay of their distance
to the cell (Canter et al. 2000; Levine 2009); `mode = "lik"` gives the
sum of the logarithms, the log-likelihood of the cell as base when the
decay function is the planar density of an incident about the base and
the incidents are independent (O’Leary 2009), so that `prob` is the
posterior over cells under a uniform prior. The sum surface is not a
likelihood.

``` r

s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
            ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
s$peak
#>  cell     x     y 
#> 467.0   3.3   3.1
```

## Empirical distance decay from solved trips

[`jtcal()`](https://rootcoder007.github.io/crimlab/reference/jtcal.md):
A kernel density estimate of trip lengths from known offender trips, to
calibrate a decay function (Levine 2009; Silverman 1986). With
`refl = TRUE` the kernels are reflected at 0, so no mass falls below 0
and the estimate integrates to 1 over \[0, Inf). Tabulated at evenly
spaced distances it is a `"table"` decay for
[`jtdecay()`](https://rootcoder007.github.io/crimlab/reference/jtdecay.md)
and
[`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md).

``` r

set.seed(1)
trips <- rexp(200, 1 / 2)
jtcal(trips, t = seq(0, 10, 0.5), h = 0.5)
#>  [1] 0.341941995 0.334326174 0.303888364 0.263429537 0.223093541 0.178452942 0.129191709
#>  [8] 0.095269714 0.079083025 0.063478487 0.046114240 0.036675173 0.025708978 0.013938566
#> [15] 0.006084284 0.003949455 0.004978870 0.005423329 0.005837164 0.005507356 0.003471230
```

## Accuracy of a journey-to-crime surface on a solved series

[`jtacc()`](https://rootcoder007.github.io/crimlab/reference/jtacc.md):
Given a surface over the grid of
[`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md)
and the offender’s known base: the surface’s value in the base’s cell,
the hit score (the share of the grid searched, in decreasing order of
the surface, up to and including the base’s cell), in the best and the
worst order of ties (Rossmo 2000; Canter et al. 2000), and the error
distance from the peak cell’s centre to the base.

``` r

s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
            ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
jtacc(s$prob, win = c(0, 6, 0, 6), base = c(3.2, 2.9))
#>            p      hitbest     hitworst          err         peak         cell 
#>   0.00932350   0.01000000   0.01111111   0.22360680 467.00000000 437.00000000
```

## Bayesian journey to crime

[`bjtc()`](https://rootcoder007.github.io/crimlab/reference/bjtc.md):
Combines a journey-to-crime surface over origin zones with an
origin-destination matrix of solved offences (Levine 2009; Levine and
Block 2011): `P(Jtc)`, the surface; `P(O)`, the origins of all solved
offences; `P(O|Jtc)`, the origins of offences in the zones where the
series happened; their product, Bayes’ numerator with the surface as
prior; and the risk, the product over `P(O)`, renormalised. Each block
sums to 1. The matrix describes past solved offences and inherits their
selection.

``` r

od <- matrix(c(5, 1, 0, 2, 4, 1, 0, 1, 6), 3)
bjtc(pj = c(0.5, 0.3, 0.2), od = od, dest = c(1, 2, 2))
#>    pj   po  pod       prod       risk
#> 1 0.5 0.35 0.45 0.59210526 0.55900621
#> 2 0.3 0.30 0.45 0.35526316 0.39130435
#> 3 0.2 0.35 0.10 0.05263158 0.04968944
```

## Circle hypothesis of a crime series

[`crchyp()`](https://rootcoder007.github.io/crimlab/reference/crchyp.md):
Canter and Larkin’s (1993) circle hypothesis: draw the circle whose
diameter joins the two crime sites of a series that lie farthest apart.
An offender whose home lies within it is a “marauder” (who travels out
from home), one whose home lies outside a “commuter”. Gives the circle,
the share of the sites it holds (all of them, unless the series is far
from round) and, for each home given, whether it lies within.

``` r

crchyp(c(0, 4, 2), c(0, 0, 1), home = c(2, -1))
#> $centre
#> x y 
#> 2 0 
#> 
#> $r
#> [1] 2
#> 
#> $pair
#> [1] 1 2
#> 
#> $share
#> [1] 1
#> 
#> $mar
#> [1] TRUE
```

## Hit score of a geographic profile

[`hitscr()`](https://rootcoder007.github.io/crimlab/reference/hitscr.md):
A geographic profile ranks the cells of a search area; searched from the
highest score down, the hit score is the share of the area searched
before the offender’s home is reached (Rossmo 2000), the smaller the
better. Cells scored the same as the home cell are counted both ways:
`hit` searches them after the home, `hitt` before it.

``` r

hitscr(c(0.9, 0.5, 0.7, 0.5, 0.1), home = 2)
#>   home above ties rank hit hitt
#> 1    2     2    2    3  40   80
```

## Spatial weights

[`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md):
Builds the neighbours and weights the autocorrelation statistics
([`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md),
[`gorg()`](https://rootcoder007.github.io/crimlab/reference/gorg.md),
[`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md),
[`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md))
use: the `k` nearest neighbours of each point (`type = "knn"`, ties
broken by the lower row number), the points within distance `d`
(`"band"`, the distance itself included), or a contiguity list `nb`
(`"contig"`, one element per unit giving its neighbours’ numbers, as
spdep’s `poly2nb()` returns). Weights are 1, or the inverse distance to
the power `idw`; `style = "row"` divides each unit’s weights by their
sum (a unit without neighbours keeps none), `"raw"` leaves them.

``` r

set.seed(1)
p <- cbind(x = runif(20), y = runif(20))
w <- spwt(p, k = 3)
str(w)
#> List of 5
#>  $ n: int 20
#>  $ i: int [1:60] 1 1 1 2 2 2 3 3 3 4 ...
#>  $ j: int [1:60] 12 16 19 5 11 14 13 16 17 6 ...
#>  $ d: num [1:60] 0.347 0.353 0.24 0.179 0.317 ...
#>  $ w: num [1:60] 0.333 0.333 0.333 0.333 0.333 ...
spwt(nb = list(2, c(1, 3), 2), type = "contig")$w
#> [1] 1.0 0.5 0.5 1.0
```

## Geary’s C

[`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md):
Geary’s (1954) contiguity ratio of values on units with spatial weights
`w`: C = (n - 1) sum w_ij (x_i - x_j)^2 / (2 S0 sum (x_i - xbar)^2), 1
under no autocorrelation, below 1 when neighbours are alike. Its
variance is given under normality and under randomisation (Cliff and Ord
1981), with z statistics `(1 - C) / sd` and their upper normal tails,
and, with `nsim > 0`, a permutation p-value against positive
autocorrelation. That C has expectation 1 under randomisation is proved
in Lean: see `proof("gearc")`.

``` r

set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
gearc(p[, 1] + p[, 2] + rnorm(50, sd = 0.1), spwt(p, k = 4))
#>            c           vn           vr           zn           zr           pn 
#> 9.221870e-02 9.592157e-03 9.036261e-03 9.268791e+00 9.549637e+00 9.413378e-21 
#>           pr        pperm           s0 
#> 6.507588e-22          NaN 5.000000e+01
```

## The general G of Getis and Ord

[`gorg()`](https://rootcoder007.github.io/crimlab/reference/gorg.md):
The global concentration statistic of Getis and Ord (1992) for
non-negative values: G = sum_ij w_ij x_i x_j / sum\_{i != j} x_i x_j,
its expectation `S0 / (n (n - 1))` and variance under randomisation, the
z statistic, its upper normal tail (high values near high values) and,
with `nsim > 0`, a permutation p-value. Unlike Moran’s I and Geary’s C
it tells high-value clusters from low-value ones. Its expectation under
randomisation is proved in Lean: see `proof("gorg")`.

``` r

set.seed(1)
p <- cbind(x = runif(50), y = runif(50))
gorg(exp(3 * p[, 1]), spwt(p, type = "band", d = 0.2, style = "raw"))
#>            g            e            v            z            p        pperm 
#> 1.663958e-01 1.167347e-01 1.514881e-04 4.034852e+00 2.731842e-05          NaN
```

## Local Moran (LISA)

[`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md):
Anselin’s (1995) local indicator of spatial association for each unit:
with `z = x - mean(x)` and `m2 = mean(z^2)`, I_i = z_i / m2 sum_j w_ij
z_j, its expectation and variance under total randomisation, the
two-sided normal p-value, and a conditional permutation p-value: unit
`i` keeps its value and its neighbours’ values are drawn without
replacement from the others, `nsim` times from `seed`. The quadrant
tells high-high and low-low clusters from high-low and low-high
outliers. The p-values are per unit and not adjusted for the many tests
made (use [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html));
with row-standardised weights the sum of the `I_i` is `S0` times global
Moran’s I. That identity is proved in Lean: see `proof("lmorn")`.

``` r

set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- ifelse(p[, 1] > 0.7 & p[, 2] > 0.7, 5, 1) + rnorm(60, sd = 0.3)
l <- lmorn(v, spwt(p, k = 5), nsim = 199, seed = 1)
head(l)
#>            ii           e         v          z         p pperm quad
#> 1 -0.08136911 -0.01694915 0.1668323 -0.1577177 0.8746793 0.290    4
#> 2  0.12229676 -0.01694915 0.1668323  0.3409120 0.7331698 0.215    2
#> 3  0.03723394 -0.01694915 0.1668323  0.1326550 0.8944662 0.480    2
#> 4  0.14258312 -0.01694915 0.1668323  0.3905786 0.6961088 0.095    2
#> 5  0.18202573 -0.01694915 0.1668323  0.4871449 0.6261557 0.255    2
#> 6  0.09264972 -0.01694915 0.1668323  0.2683280 0.7884469 0.390    2
```

## Local G and G\* of Getis and Ord

[`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md):
The local statistics of Getis and Ord (1992) as z values (Ord and Getis
1995): `G_i` compares the weighted sum of the neighbours’ values with
its expectation from the other `n - 1` values; G\*\_i (`star = TRUE`)
includes the unit itself with weight 1 and uses all `n`. Large positive
values mark concentrations of high values, large negative ones of low
values. With `star = TRUE` the weights are best binary
(`style = "raw"`).

``` r

set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- rpois(60, ifelse(p[, 1] > 0.7, 6, 2))
locgi(v, spwt(p, type = "band", d = 0.2, style = "raw"))
#>  [1]  0.23994257 -1.19393804 -0.85693776  4.16902467 -1.51580855  2.62650766  3.78952137
#>  [8]  0.87555456 -0.91284462 -1.74706652 -2.08423675 -0.89182354  1.19393804 -1.51232152
#> [15]  1.27896346 -0.37705262  2.28529782  2.46318889 -0.23878761  2.08423675  4.04305111
#> [22] -1.05904046  0.55717109 -2.22803818 -0.07959587 -0.72460663 -1.06343180 -0.56842821
#> [29]  2.08423675 -0.73778305 -0.72460663 -0.91284462 -0.72460663 -0.75790427  2.30828021
#> [36]  1.28199634  1.78732733 -1.97494476  2.30828021 -0.56071512  2.09530386  1.19393804
#> [43]  1.03474630  0.30383766 -1.72790811  1.03289627 -1.17365736 -1.62693646  1.44921325
#> [50]  0.78034560 -0.65203187  2.50225827 -1.44605478 -2.08423675 -1.51918828 -1.17365736
#> [57] -1.03474630 -0.91485098 -0.23878761 -0.55738971
```

## Directional mean and circular variance

[`dirmn()`](https://rootcoder007.github.io/crimlab/reference/dirmn.md):
The mean direction of a set of directions, the mean resultant length `R`
(from 0, directions spread evenly, to 1, all the same), the circular
variance `1 - R` and the circular standard deviation sqrt(-2 log R)
(Mardia and Jupp 2000). Directions are given as angles, or as points
seen from an origin (the directions of incidents from an offender’s
base, say). The mean direction of points depends on the origin chosen.
That `R` lies in \[0, 1\] is proved in Lean: see `proof("dirmn")`.

``` r

dirmn(c(1, 2, 0), c(1, 1, 3), origin = c(0, 0))
#>      theta    bearing          r        var         sd          w          n 
#>  0.9315290 36.6273216  0.8948002  0.1051998  0.4714971  3.0000000  3.0000000
dirmn(angle = c(0.1, 0.2, 6.2))
#>        theta      bearing            r          var           sd            w 
#>  0.072363804 85.853859466  0.993136769  0.006863231  0.117361758  3.000000000 
#>            n 
#>  3.000000000
```

## Local Getis-Ord Gi\* on a grid of counts

[`gstr()`](https://rootcoder007.github.io/crimlab/reference/gstr.md):
The local statistic Gi\* of Getis and Ord (1992; Ord and Getis 1995) for
each cell of a grid of counts, with binary weights on the cells within
`rad` cell widths (`rad = 1.5`: the 3 by 3 block, cell included), its
normal p-value and Benjamini-Hochberg q-value, a conditional permutation
p-value, and an exposure-aware version: given the total, the count of a
neighbourhood is binomial with the neighbourhood’s share of the exposure
under complete spatial randomness, so cells overhanging the study area
neither inflate the mean nor get flagged. No threshold is applied.

``` r

set.seed(1)
m <- matrix(rpois(100, 2), 10)
m[4:5, 4:5] <- 9
g <- gstr(m)
g[g$qe < 0.05, ]
#>    row col       gz           pz           qz  pp       ze           pe           qe
#> 34   4   4 4.374058 1.219581e-05 0.0006097906 NaN 5.368464 2.971884e-06 0.0001485942
#> 35   5   4 4.186330 2.835012e-05 0.0009450040 NaN 5.138058 7.201663e-06 0.0002400554
#> 44   4   5 4.374058 1.219581e-05 0.0006097906 NaN 5.368464 2.971884e-06 0.0001485942
#> 45   5   5 3.998602 6.371769e-05 0.0015929423 NaN 4.907652 1.698365e-05 0.0004245913
```

## Andresen’s spatial point pattern test

[`sptst()`](https://rootcoder007.github.io/crimlab/reference/sptst.md):
Compares where two point patterns fall over the same areal units
(Andresen 2009): for each unit, the difference between its share of the
test pattern `b` and its share of the base pattern `a`, with a bootstrap
p-value from multinomial resamples of each pattern; the S index is the
share of units whose difference is not significant (1: the patterns are
alike everywhere).

``` r

set.seed(1)
a <- rpois(30, 10)
b <- rpois(30, 10)
sptst(a, b, nboot = 199, seed = 1)$glob
#>        left       empty        elig     flagged           s           z           p 
#>   0.0000000   0.0000000  30.0000000   2.0000000   0.9333333  -0.4188539   0.4464579 
#>       nboot 
#> 199.0000000
```
