# Fairness of a lineup from mock witnesses

Asks whether a lineup singles out its suspect. People who never saw the
crime (mock witnesses) read the witness's description and pick the
lineup member who fits it best; in a fair lineup of `k` they pick each
member about equally often. A choice of 0 records a mock witness who
chose no one. With `n` mock witnesses, `m` of whom chose a member, and
`o_i` the choices of member `i`:

## Usage

``` r
lnsize(choice, k = NULL, susp = 1, level = 0.95)
```

## Arguments

- choice:

  The lineup position each mock witness picked: whole numbers from 0, 0
  for a mock witness who chose no one; one member at least must be
  chosen.

- k:

  `NULL` for the largest position picked, or the lineup's size: one
  whole number from 2.

- susp:

  The suspect's position: one whole number from 1 to `k`.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `lnsize`, a list: `count`, the choices of each
position; `stat`, a named vector (`n`, `k`, `nm` (the mock witnesses who
chose a member), `psusp`, `pbin`, `fsize`, `esize`, `esizeb`, `etr`,
`etrlo`, `etrhi`); and `susp`.

## Details

- `psusp`, the share of the `n` who picked the suspect, and `pbin`, the
  binomial chance of that many or more among the `m` if each member were
  picked with chance `1 / k` (a small `pbin` is bias against the
  suspect);

- `fsize`, the functional size of Wells, Leippe and Ostrom (1979), `n`
  over the number who picked the suspect (`k` for a fair lineup);

- `esize`, the effective size of Malpass (1981),
  `ka - sum |o_i - ea| / (2 ea)`: the members nobody chose (null foils)
  are left out, `ka` is the number of members chosen at least once,
  `ea = m / ka`, and the sum runs over those `ka` members; it reduces
  the size by how far members are chosen below chance. `esizeb` is the
  same without the adjustment (`k` and `e = m / k`, the sum over all
  members), the form Tredoux (1996, Chapter 6, equation 5; 1999, p. S13)
  prefers because a member can go unchosen by chance;

- `etr`, the effective size `E' = 1 / (1 - I) = 1 / sum p_i^2` of
  Tredoux (1998), `p_i = o_i / m`, with its interval: the large-sample
  interval of `I = 1 - sum p_i^2`, whose variance is
  `4 (sum p^3 - (sum p^2)^2) / m` (Agresti and Agresti 1978), with both
  ends taken to `1 / (1 - I)` (Tredoux 1996, Chapter 6, equations 8, 9
  and 11). `etrhi` is `Inf` when the upper end of `I` reaches 1.

The tests reproduce Tables 4.1 (p. 67), 4.2 (p. 70) and 6.6 (p. 144) of
Tredoux (1996) and Table 1 (p. S10) of Tredoux (1999). Two printed
values there do not follow from their own counts and formulas: lineup i
of Table 4.2 (1, 0, 12, 12, 0, 11) has an adjusted effective size of
3.11, printed 3.17 (the unadjusted value), and the lower end of `E'` in
Table 6.6 is 1 / (1 - 0.7015) = 3.350, printed 3.503.

## References

Malpass, R. S. (1981). Effective size and defendant bias in eyewitness
identification lineups. *Law and Human Behavior*, 5(4), 299-309.
[doi:10.1007/BF01044945](https://doi.org/10.1007/BF01044945)

Tredoux, C. G. (1996). *Evaluating the Fairness of Identification
Parades with Measures of Facial Similarity*. PhD thesis, University of
Cape Town.

Tredoux, C. G. (1998). Statistical inference on measures of lineup
fairness. *Law and Human Behavior*, 22(2), 217-237.
[doi:10.1023/A:1025746220886](https://doi.org/10.1023/A%3A1025746220886)

Tredoux, C. (1999). Statistical considerations when determining measures
of lineup size and lineup bias. *Applied Cognitive Psychology*, 13(S1),
S9-S26.
[doi:10.1002/(SICI)1099-0720(199911)13:1+\<S9::AID-ACP634\>3.0.CO;2-1](https://doi.org/10.1002/%28SICI%291099-0720%28199911%2913%3A1%2B%3CS9%3A%3AAID-ACP634%3E3.0.CO%3B2-1)

Agresti, A. and Agresti, B. F. (1978). Statistical analysis of
qualitative variation. *Sociological Methodology*, 9, 204-237.
[doi:10.2307/270810](https://doi.org/10.2307/270810)

Wells, G. L., Leippe, M. R. and Ostrom, T. M. (1979). Guidelines for
empirically assessing the fairness of a lineup. *Law and Human
Behavior*, 3(4), 285-293.
[doi:10.1007/BF01039807](https://doi.org/10.1007/BF01039807)

## Examples

``` r
# 40 mock witnesses, a lineup of six, the suspect in position 3
ch <- c(rep(3, 15), rep(1, 5), rep(2, 6), rep(4, 4), rep(5, 7), rep(6, 3))
lnsize(ch, k = 6, susp = 3)
#> lnsize: 40 mock witnesses, lineup of 6, suspect in position 3
#>            n            k           nm        psusp         pbin        fsize 
#> 40.000000000  6.000000000 40.000000000  0.375000000  0.001282306  2.666666667 
#>        esize       esizeb          etr        etrlo        etrhi 
#>  4.700000000  4.700000000  4.444444444  3.345535088  6.618389805 
```
