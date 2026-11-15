# Atmospheric dispersion: sigmas, plume rise, plume, puff, grid and particles

Six standard models of how a pollutant released into the air spreads
downwind, for exposure estimates around a stack, a fire or an industrial
site (for example when relating air quality to health or justice
outcomes by neighbourhood):

- `pgsig()`: how wide and how tall a plume has grown at a given
  distance, by atmospheric stability class;

- `briggs()`: how far a hot plume rises above the stack top;

- `gplume()`: the steady concentration at receptor points from a
  continuous release;

- `gspuff()`: the concentration from a single instantaneous release,
  some time after it;

- `advdif2d()`: a concentration field on a grid moved by wind and spread
  by diffusion, step by step;

- `lagrprt()`: the same physics as a random walk of many particles.

They are screening models with textbook assumptions (flat terrain,
steady uniform wind, no chemistry or deposition), not regulatory model
runs.

In criminology these models serve environmental criminology in its two
established senses. In the first, pollution is itself the harm or the
offence: the study of environmental crime and harm (green criminology)
asks who bears the exposure from an emitting site, and the modelled
concentrations by area and group are the exposure behind burden and
environmental justice measures
([`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md),
[`polbrdarea`](https://rootcoder007.github.io/crimlab/reference/polbrdarea.md),
[`expconcidx`](https://rootcoder007.github.io/crimlab/reference/expconcidx.md)).
In the second, the environment is a driver of offending and
victimisation: studies relating air pollution to crime join an exposure
series to crime counts by place and day (Bondy, Roth and Sager 2020, on
London wards; Herrnstadt et al. 2019, on Los Angeles), and these models
supply the exposure where no monitor stands.

## Usage

``` r
pgsig(x, stab = "D", set = "rural")

briggs(
  x,
  u,
  diam,
  exitvel,
  stktmp,
  ambtmp,
  stab = "D",
  dthdz = NULL,
  g = 9.80616
)

gplume(q, u, h, rcpt, stab = "D", set = "rural", mixh = NULL, nimg = 3)

gspuff(mass, u, h, rcpt, t, stab = "D", set = "rural", mixh = NULL, nimg = 3)

advdif2d(c0, u, v, kx, ky, dx, dy, dt, nsteps, source = NULL, maxwork = 1e+08)

lagrprt(
  n,
  x0,
  y0,
  u,
  v,
  kx,
  ky,
  dt,
  nsteps,
  seed = 0,
  grid = NULL,
  normals = NULL,
  maxwork = 1e+08
)

# S3 method for class 'advdif2d'
print(x, ...)

# S3 method for class 'lagrprt'
print(x, ...)
```

## Arguments

- x:

  Downwind distance(s) in metres.

- stab:

  Pasquill stability class, one of `"A"` to `"F"`.

- set:

  `"rural"` (the default) or `"urban"`: which set of Briggs coefficients
  to use.

- u:

  Wind speed in m/s (the x velocity for the grid models).

- diam, exitvel, stktmp, ambtmp:

  Stack diameter (m), exit velocity (m/s), stack gas and ambient
  temperatures (K). Only buoyant rise is modelled: a stack at the air's
  temperature (within half a kelvin) has zero rise, and momentum rise
  from the exit velocity is not added.

- dthdz:

  Potential temperature gradient (K/m) for the stable classes; the
  defaults are 0.020 (E) and 0.035 (F).

- g:

  Gravitational acceleration in m/s^2, default 9.80616.

- q:

  Emission rate (mass per second, for example g/s).

- h:

  Effective stack height (m): physical height plus plume rise.

- rcpt:

  A matrix with columns x, y, z, or a list of length-3 vectors, or a
  vector of x, y, z triples.

- mixh:

  Mixing-lid height in metres, or `NULL` (the default) for no lid. Under
  a lid the plume is reflected between the ground and the lid.

- nimg:

  Whole number, default 3: the number of lid image pairs (0 to 1,000).
  The terms decay quickly; a few are enough.

- mass:

  Puff mass (for example g).

- t:

  Time since release (s).

- c0:

  Initial concentration grid: a non-negative numeric matrix, one row per
  x cell and one column per y cell.

- v:

  Wind velocity along y (m/s), may be negative.

- kx, ky:

  Eddy diffusivities along x and y (m^2/s), non-negative.

- dx, dy:

  Grid spacings (m), positive.

- dt:

  Time step (s), positive.

- nsteps:

  Whole number of time steps (0 returns the start).

- source:

  Source grid (concentration per unit time), or `NULL`.

- maxwork:

  The most cell updates (cells x steps x sub-steps; for `lagrprt()`,
  particles x steps) a run may take, 100 million by default; `Inf` lifts
  it.

- n:

  Whole number of particles, at most 100 million.

- x0, y0:

  Release point (m).

- seed:

  Number, default 0: the seed of the particle random walk, a whole
  number from -2147483647 to 2^53 - 1.

- grid:

  `NULL`, or `c(xmin, xmax, ymin, ymax, nx, ny)` to count particles into
  a concentration grid.

- normals:

  `NULL` (the default stream above), or a function `(n, seed, stream)`
  returning `n` standard normals; step `k` (from 0) draws streams `2k`
  for x and `2k + 1` for y.

- ...:

  Not used.

## Value

`pgsig()`: a list with `sigmay` and `sigmaz`. `briggs()`: a list with
`flux`, `finalrise`, `xfinal` and `rise` (at each `x`). `gplume()` and
`gspuff()`: a numeric vector of concentrations, one per receptor.
`advdif2d()`: a list of class `advdif2d` with `field`, `mass`, `cfl`,
`difn`, `stabn`, `stable`, `dt`, `substeps`, the axes `x` and `y`, `dx`,
`dy`, `units` and `title`. `lagrprt()`: a list of class `lagrprt` with
`x`, `y` (final particle positions), `mean_x`, `meany`, `grid` (as
given) and, with a grid, `concentration` (an `nx` by `ny` matrix of
particle density per unit area, summing to the share of particles inside
the grid divided by the cell area).

## Details

`pgsig()`: Briggs (1973) fits of the Pasquill-Gifford dispersion
coefficients by stability class A-F, rural
`sigmay = a_y x (1 + 0.0001 x)^(-1/2)` (urban: 0.0004) and
`sigmaz = a_z x (1 + b_z x)^(e_z)`. `briggs()`: buoyancy flux
`F = g v_s d^2 (T_s - T_a) / (4 T_s)`; final rise `21.425 F^(3/4) / u`
with `x_f = 49 F^(5/8)` when `F < 55`, else `38.71 F^(3/5) / u` with
`x_f = 119 F^(2/5)`, for classes A-D; `2.6 (F / (u s))^(1/3)` with
`s = g (dtheta/dz) / T_a` and `x_f = 2.0715 u / sqrt(s)` for E and F;
transitional rise `1.6 F^(1/3) x^(2/3) / u` before `x_f` (EPA ISC3).
`gplume()`: `Q / (2 pi u sigmay sigmaz) exp(-y^2 / (2 sigmay^2)) V`,
with `V` the vertical term with ground reflection and, under a mixing
lid, `nimg` pairs of image terms. `gspuff()`: an instantaneous release
with `sigma_x = sigmay` evaluated at the travel distance `u t`.
`advdif2d()`: explicit upwind advection, central diffusion, forward
Euler, zero boundaries. The explicit scheme is stable when
`|u dt/dx| + |v dt/dy| + 2 (kx dt/dx^2 + ky dt/dy^2) <= 1` (the von
Neumann bound for upwind advection with forward-time centred-space
diffusion); that sum is returned as `stabn` with `stable`. A call
outside the bound is not run as given: each step is split into enough
sub-steps of a smaller `dt` to satisfy it (reported as `substeps` and
`dt`), with a warning, so the field is a solution over the same physical
time rather than a numerical explosion. The separate `cfl` and `difn`
are reported for reference; each below 1 is necessary, not sufficient.
`lagrprt()`: the random walk `x + u dt + sqrt(2 K dt) xi` (Thomson 1987)
with an optional counting grid. The standard normals `xi` of step `k`
(from 0) come from the splitmix64 uniform stream (Steele, Lea and Flood
2014) seeded with `(seed * 65536 + s) mod 2147483647`, `s = 2k` for x
and `2k + 1` for y, through the normal quantile; R's own random number
stream is not touched.

The stability classes run from `"A"` (very unstable, sunny with light
wind) through `"D"` (neutral, overcast or windy) to `"F"` (stable, clear
night with light wind). Units are SI: metres, seconds, kelvin; with `q`
in g/s, `gplume()` returns g/m^3, and with `mass` in g, `gspuff()`
returns g/m^3. The plume is not defined at or upwind of the source: a
receptor with \\x \le 0\\ gets 0. The plume's `h` is the effective
height, the stack height plus the rise from `briggs()`.

Every physical input is checked: emission rates, masses, heights and
distances are non-negative, wind speed, travel time and grid steps
positive, `nimg` a whole number up to 1,000, the stack as hot as the air
or hotter (temperatures below 150 are refused as probably Celsius), and
`x0`, `y0` single numbers; a release above `mixh` warns, and so does a
receptor below ground or above the lid (the reflection terms give mirror
values there). Unknown stability classes or settings are errors. Both
grid models refuse a run above `maxwork` cell or particle updates.

## Boundary

the outermost ring of cells is held at zero (an open boundary): mass
that reaches it leaves the domain and `mass` falls accordingly. Pad `c0`
with empty cells so the plume of interest stays inside the ring.

## References

Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
contemporaneous relationship between air pollution and crime. Journal of
the Association of Environmental and Resource Economists 7(3), 555-585.
[doi:10.1086/707127](https://doi.org/10.1086/707127) . Herrnstadt, E.,
Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air pollution as a
cause of violent crime. Working paper. Briggs, G. A. (1973). Diffusion
estimation for small emissions. ATDL Contribution 79. Briggs, G. A.
(1975). Plume rise predictions. American Meteorological Society. Turner,
D. B. (1970). Workbook of Atmospheric Dispersion Estimates. EPA AP-26.
Seinfeld, J. H. and Pandis, S. N. (2016). Atmospheric Chemistry and
Physics, 3rd ed. Thomson, D. J. (1987). Criteria for the selection of
stochastic models of particle trajectories in turbulent flows. Journal
of Fluid Mechanics 180, 529-556. Steele, G. L., Lea, D. and Flood, C. H.
(2014). Fast splittable pseudorandom number generators. OOPSLA 2014,
453-472.
[doi:10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195)

## Examples

``` r
# plume widths 1 km downwind: neutral, then very stable, urban
pgsig(1000, "D")
#> $sigmay
#> [1] 76.27701
#> 
#> $sigmaz
#> [1] 37.94733
#> 
pgsig(c(500, 1000, 2000), "F", set = "urban")
#> $sigmay
#> [1]  50.20790  92.96697 163.97832
#> 
#> $sigmaz
#> [1] 30.23716 50.59644 80.00000
#> 

# rise of a hot plume, and the effective height it gives a 50 m stack
rise <- briggs(c(100, 1000), u = 4, diam = 1, exitvel = 12,
               stktmp = 420, ambtmp = 285)
rise$rise
#> [1] 18.22335 28.88281
heff <- 50 + rise$finalrise

# ground-level concentration (g/m^3) from 100 g/s along the centreline
rec <- cbind(x = c(500, 1000, 2000, 8000), y = 0, z = 0)
gplume(100, u = 4, h = heff, rcpt = rec)
#> [1] 2.120279e-05 3.168716e-04 3.826280e-04 1.051326e-04

# under a 120 m mixing lid the far receptor sees a trapped plume
gplume(100, u = 4, h = heff, rcpt = rec, mixh = 120)
#> [1] 2.120279e-05 3.172063e-04 4.073068e-04 1.738501e-04

# a 1 kg puff, 10 minutes after release, at the puff centre
gspuff(1000, u = 3, h = 10, rcpt = rbind(c(1800, 0, 0)), t = 600)
#> [1] 0.0001266794

# a small grid: a blob moved by wind and spread by diffusion
c0 <- matrix(0, 15, 15)
c0[4, 8] <- 1
f <- advdif2d(c0, u = 1, v = 0, kx = 0.5, ky = 0.5,
              dx = 1, dy = 1, dt = 0.2, nsteps = 20)
f$stable
#> [1] TRUE
round(f$mass, 3)
#> [1] 0.962

# 2,000 simulated particles, counted on a grid
p <- lagrprt(2000, x0 = 0, y0 = 0, u = 1, v = 0, kx = 1, ky = 1,
             dt = 1, nsteps = 10, seed = 1,
             grid = c(-10, 30, -20, 20, 8, 8))
c(p$mean_x, p$meany)
#> [1] 10.033977334 -0.003237509

# inputs that cannot be right are refused
try(pgsig(100, "G"))
#> Error : pgsig: `stab` must be one of the Pasquill classes "A" to "F".
try(briggs(100, u = 4, diam = 1, exitvel = 12, stktmp = 155, ambtmp = 20))
#> Error : briggs: temperatures are in kelvin: stktmp 155 and ambtmp 20 look like degrees Celsius (ambient air is about 290 K).
```
