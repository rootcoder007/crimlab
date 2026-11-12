# Co-offending graph from an incident table

The one-mode projection of offenders onto incidents: two offenders are
joined when they are recorded in the same incident, and the weight of
the tie is the number of incidents they share. An offender recorded
twice in one incident counts once there. Offenders who never offend with
anyone appear in `nodes` but in no edge; give `nodes` to the network
functions to keep them as isolated nodes.

## Usage

``` r
coofg(inc, off)
```

## Arguments

- inc:

  Incident identifiers, one per row: a vector of labels with no missing
  values.

- off:

  Offender identifiers, one per row: labels, same length.

## Value

A list: `edges`, a data frame with `from`, `to` (offender labels) and
`w` (incidents shared), and `nodes`, every offender in order of first
appearance.

## Details

Ties are listed in the order each pair is first met: incidents in the
order they first appear, and within an incident the pairs of its
offenders in the order of the rows; `from` is the offender who appears
first in `off`. Labels are compared as text
([`as.character()`](https://rdrr.io/r/base/character.html); a factor by
its labels).

## References

Reiss, A. J. (1988). Co-offending and criminal careers. *Crime and
Justice*, 10, 117-170.
[doi:10.1086/449145](https://doi.org/10.1086/449145)

## Examples

``` r
g <- coofg(inc = c(1, 1, 1, 2, 2, 3), off = c("a", "b", "c", "a", "b", "d"))
g$edges
#>   from to w
#> 1    a  b 2
#> 2    a  c 1
#> 3    b  c 1
ntdeg(g$edges$from, g$edges$to, w = g$edges$w, nodes = g$nodes)
#>   node deg wdeg
#> 1    a   2    3
#> 2    b   2    3
#> 3    c   2    2
#> 4    d   0    0
```
