# Accuracy of a hot-spot map

Asks how well a map of predicted hot spots (from
[`kdsurf()`](https://rootcoder007.github.io/crimlab/reference/kdsurf.md),
[`prohs()`](https://rootcoder007.github.io/crimlab/reference/prohs.md),
[`rtmfit()`](https://rootcoder007.github.io/crimlab/reference/rtmfit.md)
or any score per cell) caught the crimes that followed. The hot cells
are those of highest predicted score whose areas together first reach
the share `top` of the study area (or the cells marked `TRUE` when
`pred` is logical). Then:

## Usage

``` r
hsacc(pred, obs, area = NULL, top = 0.1, prev = NULL)
```

## Arguments

- pred:

  The predicted score of each cell (higher is hotter), or `TRUE` for the
  cells of a map's hot spots.

- obs:

  The crimes of each cell in the evaluation period: numbers from 0, not
  all 0.

- area:

  `NULL` for cells of equal area, or the area of each cell: positive.

- top:

  The shares of the area to mark hot: numbers above 0 and at most 1, one
  row of the result each (not read when `pred` is logical).

- prev:

  `NULL` or the crimes of each cell in the period the map was made from:
  numbers from 0.

## Value

A data frame with one row per `top`: `top`, `cells`, `area` (hot),
`ashare`, `n`, `N`, `hit`, `pai`, `nmax`, `pei` and `rri` (`NA` without
`prev`).

## Details

- the hit rate `n / N`, the share of the crimes of the evaluation period
  (`obs`) in hot cells, and the predictive accuracy index
  `PAI = (n / N) / (a / A)`, the hit rate over the hot share of the area
  (Chainey, Tompson and Uhlig 2008);

- the predictive efficiency index `PEI = n / nmax`, the crimes caught
  over the most an area of the same size could have caught, its cells
  taken by observed crimes per unit area (Drawve 2016): 1 for a perfect
  map of that size;

- with `prev`, the counts of the period the map was made from, the
  recapture rate index of Levine (2008), `(n / n_prev) / (N / N_prev)`:
  whether the hot cells kept their crimes better than the whole area.

## References

Chainey, S., Tompson, L. and Uhlig, S. (2008). The utility of hotspot
mapping for predicting spatial patterns of crime. *Security Journal*,
21(1-2), 4-28.
[doi:10.1057/palgrave.sj.8350066](https://doi.org/10.1057/palgrave.sj.8350066)

Levine, N. (2008). The "hottest" part of a hotspot: comments on "The
utility of hotspot mapping for predicting spatial patterns of crime".
*Security Journal*, 21(4), 295-302.
[doi:10.1057/sj.2008.5](https://doi.org/10.1057/sj.2008.5)

Drawve, G. (2016). A metric comparison of predictive hot spot techniques
and RTM. *Justice Quarterly*, 33(3), 369-397.
[doi:10.1080/07418825.2014.904393](https://doi.org/10.1080/07418825.2014.904393)

## Examples

``` r
set.seed(12)
risk <- rgamma(400, 0.5)
before <- rpois(400, risk)
after <- rpois(400, risk)
hsacc(pred = before, obs = after, top = c(0.05, 0.1, 0.2), prev = before)
#>    top cells area ashare  n   N       hit      pai nmax       pei       rri
#> 1 0.05    20   20   0.05 36 177 0.2033898 4.067797   72 0.5000000 0.5649718
#> 2 0.10    40   40   0.10 58 177 0.3276836 3.276836  108 0.5370370 0.5567440
#> 3 0.20    80   80   0.20 79 177 0.4463277 2.231638  148 0.5337838 0.5386713
```
