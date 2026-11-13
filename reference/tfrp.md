# Glass fragments recovered under transfer and persistence

The probability of recovering 0, 1, ... fragments of glass from a
suspect under the transfer, persistence and recovery model of Curran,
Triggs, Buckleton, Walsh and Hicks (1998), with uncertain parameters
averaged over: the transfer amount and the distance deterministically
(by quadrature), the time since the offence exactly, and the loss and
recovery rates by `nsim` Monte Carlo draws (exactly when every rate is
fixed).

## Usage

``` r
tfrp(
  par,
  deffect = TRUE,
  tdist = "nb",
  ymax = 20,
  nsim = 1000,
  seed = NULL,
  hl = 0.02,
  nq = 256
)
```

## Arguments

- par:

  The 17 parameters, in order: `lam`, `lsd` (transfer amount N(lam,
  lsd^2), used as its absolute value), `d` (distance scale), `Q` (share
  of high-persistence fragments), `l0`, `u0` (first-hour loss of the
  others, uniform), `ls0`, `us0` (of the high-persistence ones), `lj`,
  `uj` and `lsj`, `usj` (each later hour), `lR`, `uR` (share recovered),
  `lt`, `ut` (hours since the offence: fixed when equal) and `r` (the
  negative binomial probability of the time).

- deffect:

  `TRUE` for the distance effect, `FALSE` without it.

- tdist:

  The law of the time when `lt < ut`: `"nb"` (negative binomial with
  size `(lt + ut) / 2` and probability `r`), `"nbtrunc"` (that law on
  `floor(lt)` to `ceiling(ut)`) or `"unif"` (uniform on those hours).

- ymax:

  The largest count reported: a whole number from 0.

- nsim:

  The Monte Carlo draws of the rates: a whole number from 1.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

- hl, nq:

  The quadrature's grid step (1e-4 to 1) and panels (1 to 100000); the
  defaults resolve the model's usual settings.

## Value

A list: `prob`, a data frame of `y` (0 to `ymax`), `p` and `se` (the
Monte Carlo standard error, 0 when exact); `cut`, the time mass left out
and the intensity mass on the grid.

## References

Curran, J. M., Triggs, C. M., Buckleton, J. S., Walsh, K. A. J. and
Hicks, T. (1998). Assessing transfer probabilities in a Bayesian
interpretation of forensic glass evidence. *Science and Justice*, 38(1),
15-21.
[doi:10.1016/S1355-0306(98)72068-4](https://doi.org/10.1016/S1355-0306%2898%2972068-4)

## Examples

``` r
par <- c(80, 50, 0.5, 0.05, 0.8, 0.9, 0.1, 0.2, 0.4, 0.6, 0.05, 0.1,
         0.5, 0.7, 1, 1, 0.5)
tfrp(par, ymax = 10, nsim = 100, seed = 1)$prob
#>     y          p           se
#> 1   0 0.10287150 0.0010164458
#> 2   1 0.06268975 0.0006369714
#> 3   2 0.05098681 0.0005238537
#> 4   3 0.04475979 0.0004628842
#> 5   4 0.04070553 0.0004220358
#> 6   5 0.03776832 0.0003907530
#> 7   6 0.03548773 0.0003643062
#> 8   7 0.03362473 0.0003402004
#> 9   8 0.03204053 0.0003169996
#> 10  9 0.03064810 0.0002938519
#> 11 10 0.02938989 0.0002702679
```
