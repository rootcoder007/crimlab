# Hawkes process: self-exciting event times

Fits a self-exciting (Hawkes) point process to event times, such as the
times of shootings, burglaries or calls for service, by maximum
likelihood: a background rate and, on top of it, a rise in the rate
after each event that dies away with the time since. The fitted
branching ratio `eta` is the expected number of events each event
triggers; it says how much of the record is contagion rather than
background.

## Usage

``` r
hwkfit(
  times,
  hor,
  kernel = "exp",
  base = "const",
  per = 365.25,
  start = NULL,
  res = NULL,
  seed = NULL,
  maxit = 500
)
```

## Arguments

- times:

  Event times: numbers in one unit (days, say), sorted, in `[0, hor]`;
  at least 10 of them.

- hor:

  The end of the observation window: one positive number.

- kernel:

  The triggering kernel: one or more of `"exp"`, `"weibull"` and
  `"lomax"`.

- base:

  The baseline: one or both of `"const"` and `"cycle"`.

- per:

  The length of the cycle of the `"cycle"` baseline, in the unit of
  `times` (365.25 for a yearly cycle in days).

- start:

  `NULL`, or the starting parameters of a single kernel and baseline in
  the order of `par` (`eta` between 0 and 1).

- res:

  `NULL`, or the resolution the times were recorded to: one positive
  number.

- seed:

  `NULL` to draw the seed of `res`'s spreading from R's generator, or
  one whole number from 0.

- maxit:

  The most BFGS steps per start: one whole number from 1.

## Value

An object of class `hwkfit`, a list: `par`, the estimates (`a0` or `a0`
to `a3` on the log scale, `eta`, then `beta` or `shape` and `scale`);
`se` and `cov`, their standard errors and covariance; `loglik`, `aic`,
`bic`; `branch` (`eta`, the mean cluster `size`, the stationary `rate`
of a constant baseline and the triggered `share`); `conv` (`TRUE` when
the gradient vanished) and `iter`; `kernel`, `base`, `per`, `n`, `hor`,
`jit` (whether `res` spread the times) and `times` (the times fitted);
and `cmp`, one row per kernel and baseline fitted (`kernel`, `base`,
`loglik`, `k` parameters, `aic`, `bic`, `eta`, `conv`). A branching
ratio above 0.99 is warned about.

## Details

The intensity is \$\$\lambda(t) = \nu(t) + \eta \sum\_{t_j \< t} g(t -
t_j)\$\$ over the window `[0, hor]`, with a normalised kernel `g`:
exponential (`"exp"`, rate `beta`: a Markovian process), Weibull
(`"weibull"`, `shape` and `scale`) or Lomax (`"lomax"`, the power-law
kernel of seismology's modified Omori law, `shape` and `scale`, density
`shape / scale (1 + t / scale)^-(shape + 1)`). The baseline is constant
(`"const"`, \\\nu = e^{a_0}\\) or a log-linear trend with one cycle of
length `per` (`"cycle"`, \\\nu(t) = \exp(a_0 + a_1 t / hor + a_2 \sin(2
\pi t / per) + a_3 \cos(2 \pi t / per))\\), whose integral is taken by
Simpson's rule on steps of at most `per / 16`. The log-likelihood is
\\\sum_i \log \lambda(t_i) - \int_0^{hor} \lambda\\; it is maximised
with its analytic gradient by BFGS over the baseline terms, the logit of
`eta` (so `0 < eta < 1`) and the logarithms of the kernel parameters.
The exponential kernel is fitted from two starts and the better optimum
kept; the others from their own start and from the exponential fit (the
shape-1 member). Standard errors come from the inverse of the numerical
Hessian of the analytic gradient.

Several kernels or baselines fit every combination: `cmp` compares them
(the exponential kernel against the others is the comparison of a
Markovian with a non-Markovian process) and the one with the least AIC
is returned. `branch` holds the arithmetic of the branching process: a
cluster started by one background event has `1 / (1 - eta)` events on
average, a constant baseline \\\nu\\ gives the stationary rate \\\nu /
(1 - eta)\\, and the share of events that are triggered is `eta`
(`proof("hwkfit")`).

Events at one instant do not excite each other: the sum runs over
strictly earlier events. Times recorded to a resolution (dates, whole
hours) carry ties the process never makes; `res` spreads each time
uniformly over `[t, t + res)` first (cut at `hor`), from `seed`. The
exponential kernel costs O(n) per likelihood, the others O(n^2).

## References

Hawkes, A. G. (1971). Spectra of some self-exciting and mutually
exciting point processes. *Biometrika*, 58(1), 83-90.
[doi:10.1093/biomet/58.1.83](https://doi.org/10.1093/biomet/58.1.83)

Ozaki, T. (1979). Maximum likelihood estimation of Hawkes' self-exciting
point processes. *Annals of the Institute of Statistical Mathematics*,
31, 145-155.
[doi:10.1007/BF02480272](https://doi.org/10.1007/BF02480272)

Mohler, G. O., Short, M. B., Brantingham, P. J., Schoenberg, F. P. and
Tita, G. E. (2011). Self-exciting point process modeling of crime.
*Journal of the American Statistical Association*, 106(493), 100-108.
[doi:10.1198/jasa.2011.ap09546](https://doi.org/10.1198/jasa.2011.ap09546)

## See also

[`hwkres()`](https://rootcoder007.github.io/crimlab/reference/hwkres.md)
for the residuals,
[`hwksim()`](https://rootcoder007.github.io/crimlab/reference/hwksim.md)
to simulate.

## Examples

``` r
x <- hwksim(c(a0 = log(0.5), eta = 0.5, beta = 2), 200, seed = 1)
fit <- hwkfit(x$t, 200)
fit
#> hwkfit: 219 events in [0, 200]; exp kernel, const baseline
#>             est         se
#> a0   -0.3948088 0.13638108
#> eta   0.3874376 0.07837729
#> beta  2.4507950 0.80140193
#> branching ratio 0.3874376, mean cluster size 1.632487; log-likelihood -177.8257, AIC 361.6514
hwkfit(x$t, 200, kernel = c("exp", "lomax"))$cmp
#>   kernel  base    loglik k      aic      bic       eta conv
#> 1    exp const -177.8257 3 361.6514 371.8186 0.3874376 TRUE
#> 2  lomax const -177.8102 4 363.6204 377.1767 0.3948362 TRUE
```
