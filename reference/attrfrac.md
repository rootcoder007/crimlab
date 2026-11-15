# Population attributable fraction

The share of all cases in a population that is due to an exposure: the
fraction that would not occur if the exposure were removed and nothing
else changed. It combines how strongly the exposure raises risk (the
relative risk) with how common it is (the prevalence); multiplied by the
number of cases it gives the attributable cases, as
[`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
does.

## Usage

``` r
attrfrac(rr, expoprev, rrlo = NULL, rrup = NULL)

# S3 method for class 'attrfrac'
print(x, ...)
```

## Arguments

- rr:

  Relative risk among the exposed, a single finite positive number.

- expoprev:

  Proportion of the population exposed, in `[0, 1]`.

- rrlo, rrup:

  Optional confidence limits of `rr`, both or neither.

- x:

  An `attrfrac` result.

- ...:

  Not used.

## Value

A single number of class `attrfrac`, the attributable fraction: in
`[0, 1)` for `rr >= 1`, negative for a protective exposure; when the
limits were given it carries an attribute `ci`, a named numeric vector
(`lower`, `upper`), and it keeps its inputs as attributes `rr` and
`expoprev`. [`as.numeric()`](https://rdrr.io/r/base/numeric.html) gives
the plain number.

## Details

Levin's formula (Levin 1953; Rothman, Greenland and Lash 2008, chapter
5): \$\$PAF = p (RR - 1) / (1 + p (RR - 1)),\$\$ with \\p\\ the
proportion exposed. With everyone exposed it is \\1 - 1/RR\\. It assumes
the relative risk is causal and unconfounded and that removing the
exposure changes nothing else. A relative risk below 1 gives a negative
fraction (the exposure protects).

When both confidence limits of the RR are given, the fraction's limits
are the formula applied to them: the formula increases with the RR, so
the interval of the fraction is the image of the interval of the risk.
Refusals: an `rr` or limit that is not a single finite positive number
(logical values included), a prevalence outside \\\[0, 1\]\\, only one
of the two limits, or limits that do not enclose `rr`.

## References

Rothman, K. J., Greenland, S. and Lash, T. L. (2008). Modern
Epidemiology, 3rd ed., chapter 5. Levin, M. L. (1953). The occurrence of
lung cancer in man. Acta Unio Internationalis Contra Cancrum 9, 531-541.

## See also

[`crfpm25`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md)
for a relative risk from a concentration,
[`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
for the attributable cases.

## Examples

``` r
attrfrac(1.5, 0.4)
#> Attributable fraction 0.1667 (relative risk 1.5, 0.4 exposed)
attrfrac(2, 1) # everyone exposed: 1 - 1/RR
#> Attributable fraction 0.5 (relative risk 2, 1 exposed)
paf <- attrfrac(1.08, 1, rrlo = 1.06, rrup = 1.09)
attr(paf, "ci")
#>      lower      upper 
#> 0.05660377 0.08256881 
attrfrac(0.8, 0.5) # protective
#> Attributable fraction -0.1111 (relative risk 0.8, 0.5 exposed)
#> Negative: the exposure protects.
try(attrfrac(1.5, 40))
#> Error : attrfrac: expoprev must be in [0,1]; got 40.
```
