# Tests of discrimination in searches: hit rates, thresholds, rank order

Asks whether the police search people of one group on less evidence than
people of another. Three tests read the same counts of stops, searches
and successful searches ("hits", contraband found):

## Usage

``` r
thrtst(
  stops,
  srch,
  hits,
  grp,
  area = NULL,
  off = NULL,
  ref = NULL,
  level = 0.95,
  maxit = 2000
)
```

## Arguments

- stops:

  The stops of each row: numbers from 0.

- srch:

  The searches among them: numbers from 0, at most `stops`.

- hits:

  The successful searches: numbers from 0, at most `srch`.

- grp:

  The group of the people stopped: labels, two groups at least, each
  with searches.

- area:

  `NULL` or the area (department, precinct) of each row: labels.

- off:

  `NULL` or the group of the officers of each row: labels, two groups at
  least.

- ref:

  The reference group: a label or a position; the first group by
  default.

- level:

  The confidence level: one number strictly between 0 and 1.

- maxit:

  The quasi-Newton steps of each of the three starts of each fit: one
  whole number from 1.

## Value

An object of class `thrtst`, a list: `hit`, one row per group (`group`,
`stops`, `srch`, `hits`, `srate` and `hrate`, the search and hit rates,
and `diff`, `se`, `z`, `p` of the hit rate against the reference);
`kpt`, the chi-square of equal hit rates (`chisq`, `df`, `p`); `thr`,
one row per group (`group`, `thr`, `se`, `lo`, `hi`, and `diff`, `dse`,
`z`, `p` against the reference), `NULL` when the model is not
identified; `lrt`, the likelihood-ratio test of equal thresholds
(`stat`, `df`, `p`, the log-likelihoods `llfull` and `llsame`, and the
number of `cells`); `fit`, one row per model (`model`, `conv`, `iter`);
`cells`, the fitted `phi`, `delta`, `thr`, `srate` and `hrate` of each
group and area with stops; `rank`, `NULL` without `off`, otherwise one
row per pair of officer groups and measure (`off1`, `off2`, `meas`,
`zmin`, `zmax`, `rej`); and `ref`.

## Details

- The outcome (hit-rate) test of Knowles, Persico and Todd (2001): if
  officers search whenever the chance of finding something is high
  enough, a group searched on less evidence shows a lower hit rate. The
  hit rate of each group is compared with the reference group's
  (difference, Wald standard error, z and two-sided p) and all groups
  together by Pearson's chi-square of hits and misses among the
  searched. A hit rate can differ between groups that face the same
  threshold when their risks are spread differently (the problem of
  infra-marginality), which the threshold test addresses.

- The threshold test of Simoiu, Corbett-Davies and Goel (2017), in the
  fast form of Pierson, Corbett-Davies and Goel (2017). Each stopped
  person carries a signal, normal with mean `delta` for those who carry
  contraband (a share `phi` of the stopped) and mean 0 for the rest; the
  officer's estimate of the chance of a hit is the posterior given the
  signal (a "discriminant distribution"), and a search follows when that
  chance exceeds the threshold `t`. In each group and area `phi`,
  `delta` and `t` are the logistic (`phi`, `t`) or exponential (`delta`)
  of a group term plus an area term, so a group's threshold may differ
  from another's by the same amount on the logit scale in every area.
  The model is fitted by maximum likelihood to the binomial counts of
  searches among stops and hits among searches (where the published
  method uses Bayesian inference). A lower threshold for a group is
  evidence of discrimination against it. `thr` is each group's threshold
  averaged over the areas, weighted by the areas' stops, with a
  delta-method standard error; the likelihood-ratio test compares the
  fit with one threshold for every group. The model is fitted only when
  it has more rates (two per group and area with stops) than parameters
  (`3 (groups + areas - 1)`): with `area` `NULL` it is not.

- The rank-order test of Anwar and Fang (2006), when `off` gives the
  officer's group: if officers of no group are prejudiced relative to
  the others, officers of one group search people of every group more
  often than officers of another (or less often, in every group), and
  likewise for hit rates. For each pair of officer groups and each
  measure the z statistics of the difference across the suspects' groups
  are given; `rej` is `TRUE` when one is significantly positive and
  another significantly negative at `level`, a ranking that changes with
  the suspect's group. This is a pairwise simplification of their test,
  not its joint form.

Rows are cells of a table (or single stops, with counts of 0 and 1):
rows of the same group and area (and officer group) add up.

## References

Knowles, J., Persico, N. and Todd, P. (2001). Racial bias in motor
vehicle searches: theory and evidence. *Journal of Political Economy*,
109(1), 203-229. [doi:10.1086/318603](https://doi.org/10.1086/318603)

Anwar, S. and Fang, H. (2006). An alternative test of racial prejudice
in motor vehicle searches: theory and evidence. *American Economic
Review*, 96(1), 127-151.
[doi:10.1257/000282806776157579](https://doi.org/10.1257/000282806776157579)

Simoiu, C., Corbett-Davies, S. and Goel, S. (2017). The problem of
infra-marginality in outcome tests for discrimination. *The Annals of
Applied Statistics*, 11(3), 1193-1216.
[doi:10.1214/17-AOAS1058](https://doi.org/10.1214/17-AOAS1058)

Pierson, E., Corbett-Davies, S. and Goel, S. (2017). Fast threshold
tests for detecting discrimination. arXiv:1702.08536.
[doi:10.48550/arXiv.1702.08536](https://doi.org/10.48550/arXiv.1702.08536)

## Examples

``` r
# stops, searches and hits of two groups in four areas
d <- data.frame(
  grp = rep(c("a", "b"), 4), area = rep(1:4, each = 2),
  stops = c(1200, 800, 900, 1100, 1500, 700, 1000, 1000),
  srch = c(60, 70, 50, 95, 80, 66, 52, 90),
  hits = c(24, 18, 19, 25, 30, 17, 21, 24)
)
thrtst(d$stops, d$srch, d$hits, d$grp, area = d$area)
#> thrtst: reference a
#> hit rates (chi-square 10.25 on 1 df, p 0.00137):
#>  group stops srch hits      srate     hrate       diff         se        z
#>      a  4600  242   94 0.05260870 0.3884298  0.0000000         NA       NA
#>      b  3600  321   84 0.08916667 0.2616822 -0.1267475 0.03979326 -3.18515
#>           p
#>          NA
#>  0.00144679
#> thresholds (likelihood ratio 0.03214 on 1 df, p 0.858):
#>  group          thr         se         lo        hi      diff        dse
#>      a 0.0002588271 0.08729044 -0.1708273 0.1713450 0.0000000         NA
#>      b 0.2601602016 0.02633179  0.2085508 0.3117696 0.2599014 0.09127834
#>        z           p
#>       NA          NA
#>  2.84735 0.004408487
# the rank-order test needs the officers' groups
thrtst(c(100, 100, 100, 100), c(10, 20, 30, 10), c(5, 5, 5, 5),
       grp = c(1, 1, 2, 2), off = c("x", "y", "x", "y"))$rank
#>   off1 off2   meas      zmin     zmax   rej
#> 1    x    y search -2.000000 3.651484  TRUE
#> 2    x    y    hit -1.936492 1.348400 FALSE
```
