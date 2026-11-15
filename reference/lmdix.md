# LMDI decomposition of a count over the factors of a Kaya identity

How much of a change in recorded environmental offences (pollution
cases, waste violations) between two years, or of a region's gap to a
reference region, is owed to each factor of an identity that multiplies
out to the count: the additive log-mean Divisia index (LMDI)
decomposition, temporal and spatial, as Duan et al. (2026, eqs 1-5, pp.
8-9) apply it to crimes of polluting the environment in 30 Chinese
provinces, 2016-2023.

## Usage

``` r
lmdix(
  y,
  q,
  unit,
  time,
  type = "time",
  at = NULL,
  ref = "mean",
  zero = 0,
  fnm = NULL
)

# S3 method for class 'lmdix'
plot(x, ..., main = NULL, col = NULL, pal = NULL, units = NULL)
```

## Arguments

- y:

  The count of each row (a unit in a period): numbers from 0.

- q:

  The quantities of the chain, in order: a numeric matrix or data frame
  with one row per row of `y` and one column per quantity, numbers from
  0; the column names name the factors.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

- type:

  `"time"` (change between periods) or `"space"` (gap to the reference
  region of each period).

- at:

  `NULL` or the periods to use, increasing: for `"time"`, the end points
  of consecutive intervals (default the first and last period; Duan et
  al. used 2016, 2018, 2020 and 2023); for `"space"`, the periods to
  decompose (default all).

- ref:

  `"mean"` or `"gmean"`: the reference region of `"space"`.

- zero:

  `0` (rows with a zero get `NA`) or a positive number put in place of
  the zeros.

- fnm:

  `NULL` or the names of the `ncol(q) + 1` factors.

- x:

  An `lmdix` result.

- ...:

  Passed to the base-graphics call.

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

- units:

  `NULL` (every unit) or the units whose rows are summed.

## Value

An object of class `lmdix`, a list: `dec`, one row per unit and interval
(`unit`, `from`, `to`, `y0`, `y1`, `dy`) or per unit and period (`unit`,
`time`, `ref` the reference count, `y`, `dy`), followed by one column of
contributions per factor; `fac`, the count and the factor values of
every row (`unit`, `time`, `y`, one column per factor); `ref`, for
`"space"`, the reference region's count and factors by period, else
`NULL`; `factors`, the factor names; `type`; and `refkind`.

## Details

The identity is built as a chain from the quantities `q_1, ..., q_m`
(the columns of `q`, in order):
`y = (y / q_1) (q_1 / q_2) ... (q_{m-1} / q_m) q_m`, `m + 1` factors
that multiply to `y` whatever the quantities. Duan et al.'s eq 1 is the
chain of population `P`, area `S`, emissions `W`, GDP, population again
and the unemployed `UN`:
`CPEs = CPER x PD x EIPA^-1 x PEI x EDL x UR^-1 x UN`, the case rate per
head, population density, area per unit of emission, emission per unit
of GDP, GDP per head, people per unemployed person and the unemployed
(see the example). The factor values themselves, the rates of the
identity, are returned in `fac`.

With `type = "time"` the change of each unit between consecutive periods
of `at` is `y_T - y_O = sum_j L(y_T, y_O) ln(f_jT / f_jO)` (eqs 2-3),
`L(a, b) = (a - b) / (ln a - ln b)` the logarithmic mean; a unit missing
either period is left out. With `type = "space"` each unit of each
period in `at` is set against a reference region of the same period (eqs
4-5): with `ref = "mean"`, the region whose count and quantities are the
means over the units present, so its count is the average count Duan et
al. benchmark against; with `ref = "gmean"`, the geometric means. Duan
et al. do not say how the reference region's factors were formed; the
arithmetic means of the quantities are this function's reading of "the
national average", and `"gmean"` the alternative. The decomposition has
no residual: the contributions add up to `dy`.

A zero count or quantity has no logarithm. With `zero = 0` (Duan et
al.'s "Null: no numerical value", figs 5-6) the rows it touches get `NA`
contributions; with `zero` positive, zeros are replaced by it first. The
contributions then depend on `zero`, and only slowly fade as it shrinks
(as one over its logarithm), so a replaced zero should be reported as
such.

The plot draws the contribution of each factor summed over the rows (of
`units`), as horizontal bars, and returns those sums.

## Environmental crime

The family treats offences against the environment as the outcome:
where, when and why pollution offences, illegal waste disposal and
trafficking, and safety crimes in extractive industries are recorded,
and how enforcement behaves. Its tools are `lmdix()` (decomposition over
a Kaya identity, Duan et al. 2026),
[`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md)
(two-way fixed effects with an interaction and the checks Duan et al.
report) and
[`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md)
(cross-sectional dependence of a panel, Andreatta et al. 2023), with the
package's spatial tests
([`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md),
[`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md),
[`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md))
for the spatial checks Andreatta et al. run. Recorded counts are counts
of detected, prosecuted or judged cases: Duan et al. count
first-instance judgments, Andreatta et al. the violations verified by
law enforcement and gathered by an NGO, and both say the counts carry
enforcement as well as offending (Andreatta et al. write the observed
rate as the true rate plus a measurement error that can move with
enforcement, p. 1372). Safety crimes are "violations of law by employers
that either do, or have the potential to, cause sudden death or injury
as a result of work-related activities" (Tombs and Whyte, quoted by
Pitoulis 2014, p. 2); they are mostly regulatory offences enforced by
inspectorates rather than police, mostly end in guilty pleas and fines
(p. 47), and the injuries behind them are under-reported (p. 3, note 1),
so recorded counts understate them. None of the three sources specifies
a model of recorded against true offences, and the family has none.

## References

Duan, H., Nie, C., He, B., Chen, X., Duan, Z., Chen, S., Song, J. and
Gao, Z. (2026). A multi-method framework for decoding the spatiotemporal
patterns and driving factors of environmental pollution crimes in China.
*Humanities and Social Sciences Communications*.
[doi:10.1057/s41599-026-08766-9](https://doi.org/10.1057/s41599-026-08766-9)

Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
Digging into waste: an analysis of waste crime in the Italian provinces.
*Regional Studies*, 57(7), 1367-1379.
[doi:10.1080/00343404.2022.2135697](https://doi.org/10.1080/00343404.2022.2135697)

Pitoulis, T. T. (2014). *Bloody Oil: A Critical Discourse Analysis of
Safety Crimes in the Alberta Oil and Gas Industry*. MA thesis in
criminology, University of Ottawa.

## Examples

``` r
# six provinces over four years: cases, population (millions), area
# (1000 km2), emissions (1000 t), GDP (billions) and the unemployed
# (millions), simulated
set.seed(1)
d <- expand.grid(year = c(2016, 2018, 2020, 2023), prov = LETTERS[1:6])
i <- match(d$prov, LETTERS)
d$P <- 10 * i * (1 + 0.01 * (d$year - 2016))
d$S <- 50 + 20 * i
d$W <- 30 * i * (1 - 0.03 * (d$year - 2016))
d$GDP <- 40 * i * (1 + 0.05 * (d$year - 2016))
d$UN <- 0.3 * i * (1 + 0.02 * (d$year - 2016))
d$cases <- rpois(24, 2 * i + 8 * (d$year == 2018))
q <- d[c("P", "S", "W", "GDP", "P", "UN")]
f <- c("CPER", "PD", "EIPA^-1", "PEI", "EDL", "UR^-1", "UN")
r <- lmdix(d$cases, q, d$prov, d$year, at = c(2016, 2018, 2020, 2023),
           fnm = f)
r
#> lmdix: time decomposition, 18 rows, 7 factors
#> contributions summed over the rows decomposed:
#>         dy       CPER         PD    EIPA^-1        PEI        EDL      UR^-1 
#>   1.000000  -2.821646   3.821646  13.156108 -30.253445  13.275691  -3.597145 
#>         UN 
#>   7.418791 
# the gap of each province to the average province, by year
lmdix(d$cases, q, d$prov, d$year, type = "space", fnm = f)$dec
#>    unit time       ref  y         dy       CPER         PD    EIPA^-1
#> 1     A 2016  6.833333  1 -5.8333333 -2.0307857 -2.1665160  2.1665160
#> 2     B 2016  6.833333  7  0.1666667  4.0371552 -1.8807839  1.8807839
#> 3     C 2016  6.833333  4 -2.8333333 -2.0177490 -0.3552223  0.3552223
#> 4     D 2016  6.833333 10  3.1666667  2.0561626  0.4448347 -0.4448347
#> 5     E 2016  6.833333  7  0.1666667 -2.3002157  0.9235474 -0.9235474
#> 6     F 2016  6.833333 12  5.1666667  0.2211070  1.7496734 -1.7496734
#> 7     A 2018 15.166667  8 -7.1666667  6.8691582 -7.9969646  7.9969646
#> 8     B 2018 15.166667 17  1.8333333 10.8240664 -4.3688607  4.3688607
#> 9     C 2018 15.166667 10 -5.1666667 -3.2545020 -0.8328306  0.8328306
#> 10    D 2018 15.166667 25  9.8333333  7.2060378  1.0524161 -1.0524161
#> 11    E 2018 15.166667 15 -0.1666667 -5.5464590  2.0140780 -2.0140780
#> 12    F 2018 15.166667 16  0.8333333 -7.5640268  2.9708745 -2.9708745
#> 13    A 2020  8.333333  4 -4.3333333  3.0629439 -4.2140571  4.2140571
#> 14    B 2020  8.333333  5 -3.3333333  0.3183746 -1.7744719  1.7744719
#> 15    C 2020  8.333333  8 -0.3333333  0.9253891 -0.5482282  0.5482282
#> 16    D 2020  8.333333 10  1.6666667  0.4460084  0.4889592 -0.4889592
#> 17    E 2020  8.333333  9  0.6666667 -2.4229913  1.1567012 -1.1567012
#> 18    F 2020  8.333333 14  5.6666667 -0.2206696  2.0828614 -2.0828614
#> 19    A 2023  7.000000  1 -6.0000000 -2.1372431 -2.2008204  2.2008204
#> 20    B 2023  7.000000  1 -6.0000000 -4.2744862 -0.8384777  0.8384777
#> 21    C 2023  7.000000  6 -1.0000000  0.0000000 -0.4355433  0.4355433
#> 22    D 2023  7.000000 12  5.0000000  3.7612963  0.4961877 -0.4961877
#> 23    E 2023  7.000000  8  1.0000000 -1.6710943  1.0000000 -1.0000000
#> 24    F 2023  7.000000 14  7.0000000  1.5567469  1.9257507 -1.9257507
#>              PEI           EDL         UR^-1          UN
#> 1   0.000000e+00  0.000000e+00  6.739784e-16  -3.8025477
#> 2   0.000000e+00  0.000000e+00  1.535734e-15  -3.8704886
#> 3   0.000000e+00  0.000000e+00  1.174799e-15  -0.8155843
#> 4   0.000000e+00  0.000000e+00  1.846618e-15   1.1105041
#> 5   0.000000e+00  0.000000e+00  1.535734e-15   2.4668823
#> 6   0.000000e+00  0.000000e+00  2.037369e-15   4.9455597
#> 7   0.000000e+00  2.487764e-15  0.000000e+00 -14.0358248
#> 8   0.000000e+00  3.567347e-15  0.000000e+00  -8.9907330
#> 9   0.000000e+00  0.000000e+00  0.000000e+00  -1.9121647
#> 10  0.000000e+00  4.368836e-15  0.000000e+00   2.6272955
#> 11  0.000000e+00  3.349139e-15  0.000000e+00   5.3797923
#> 12  0.000000e+00  0.000000e+00  0.000000e+00   8.3973601
#> 13 -1.310945e-15 -1.310945e-15  1.310945e-15  -7.3962773
#> 14 -1.448926e-15 -1.448926e-15  1.448926e-15  -3.6517079
#> 15  0.000000e+00 -1.813112e-15  1.813112e-15  -1.2587224
#> 16 -2.029789e-15 -2.029789e-15  2.029789e-15   1.2206583
#> 17  0.000000e+00 -1.923437e-15  1.923437e-15   3.0896580
#> 18  0.000000e+00 -2.425343e-15  2.425343e-15   5.8873362
#> 19 -3.423251e-16  0.000000e+00  0.000000e+00  -3.8627569
#> 20 -3.423251e-16  0.000000e+00  0.000000e+00  -1.7255138
#> 21 -7.202194e-16  1.440439e-15  0.000000e+00  -1.0000000
#> 22 -1.029898e-15  0.000000e+00  0.000000e+00   1.2387037
#> 23 -1.662864e-15  1.662864e-15 -8.314322e-16   2.6710943
#> 24 -1.121199e-15  2.242399e-15  0.000000e+00   5.4432531
```
