# Fairness of a lineup from mock witnesses

Asks whether a lineup singles out its suspect. People who never saw the
crime (mock witnesses) read the witness's description and pick the
lineup member who fits it best; in a fair lineup of `k` they pick each
member about equally often. From their choices:

## Usage

``` r
lnsize(choice, k = NULL, susp = 1, level = 0.95)
```

## Arguments

- choice:

  The lineup position each mock witness picked: whole numbers from 1.

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
position; `stat`, a named vector (`n`, `k`, `psusp`, `pbin`, `fsize`,
`esize`, `etr`, `etrlo`, `etrhi`); and `susp`.

## Details

- `psusp`, the share who picked the suspect, and `pbin`, the binomial
  chance of that many or more if each member were picked with chance
  `1 / k` (a small `pbin` is bias against the suspect);

- `fsize`, the functional size, the number of mock witnesses over those
  who picked the suspect (Wells, Leippe and Ostrom 1979; `k` for a fair
  lineup);

- `esize`, Malpass's (1981) effective size, `k - sum |o_i - e| / (2 e)`
  with `o_i` the choices of member `i` and `e = n / k`: how many members
  are plausible;

- `etr`, Tredoux's (1998) effective size `1 / sum p_i^2`, `p_i` the
  share choosing member `i`, with a delta-method interval from the
  multinomial variance of `sum p_i^2`, `4 (sum p^3 - (sum p^2)^2) / n`.

## References

Malpass, R. S. (1981). Effective size and defendant bias in eyewitness
identification lineups. *Law and Human Behavior*, 5(4), 299-309.
[doi:10.1007/BF01044945](https://doi.org/10.1007/BF01044945)

Tredoux, C. G. (1998). Statistical inference on measures of lineup
fairness. *Law and Human Behavior*, 22(2), 217-237.
[doi:10.1023/A:1025746220886](https://doi.org/10.1023/A%3A1025746220886)

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
#>            n            k        psusp         pbin        fsize        esize 
#> 40.000000000  6.000000000  0.375000000  0.001282306  2.666666667  4.700000000 
#>          etr        etrlo        etrhi 
#>  4.444444444  2.984575856  5.904313032 
```
