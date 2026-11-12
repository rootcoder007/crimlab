# Period-two hazards after selective depletion

A share `s` of a cohort is at high risk of reconviction, `h` a period,
and the rest at low risk, `l`. If a programme changes who survives the
first period without reconviction, the period-two hazard among the
survivors is a weighted average of `h` and `l` with weight the share of
high-risk survivors, and it rises with that share. An arm that keeps
more of its high-risk people out of reconviction in period one therefore
shows the higher period-two hazard even if it has no effect on anyone in
period two: a period-by-period hazard ratio is not a causal contrast; a
risk difference at a fixed follow-up is. See `proof("hzdsel")`.

## Usage

``` r
hzdsel(s, h, l, surh, surl)
```

## Arguments

- s:

  The high-risk share at the start, strictly between 0 and 1.

- h:

  The period hazard of the high-risk type: at most 1.

- l:

  The period hazard of the low-risk type: from 0, below `h`.

- surh:

  Period-one survival of the high-risk type in the control and treated
  arms: two numbers above 0 and at most 1, the control arm first, or
  named `control` and `treated` in either order (with `surl` named too).
  A list of single numbers or a one-row data frame is accepted.

- surl:

  Period-one survival of the low-risk type, the same way; named inputs
  are matched to `surh` by name.

## Value

A named numeric vector: `share0`, `share1` (the high-risk share among
the survivors of each arm), `haz0`, `haz1` (the period-two hazard in
each arm) and `ratio` (`haz1 / haz0`).

## References

Hernan, M. A. (2010). The hazards of hazard ratios. *Epidemiology*,
21(1), 13-15.
[doi:10.1097/EDE.0b013e3181c1ea43](https://doi.org/10.1097/EDE.0b013e3181c1ea43)

## Examples

``` r
hzdsel(s = 0.5, h = 0.5, l = 0.1, surh = c(0.5, 0.8), surl = c(0.9, 0.9))
#>    share0    share1      haz0      haz1     ratio 
#> 0.3571429 0.4705882 0.2428571 0.2882353 1.1868512 
```
