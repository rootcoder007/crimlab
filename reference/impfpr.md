# The false positive rate the other rates imply

Chouldechova's (2017) identity: in any confusion table, the false
positive rate is `p / (1 - p) (1 - ppv) / ppv (1 - fnr)`, from the base
rate `p`, the positive predictive value `ppv` and the false negative
rate `fnr`. So two groups whose scores are equally calibrated (the same
`ppv`) and miss equally (the same `fnr`) but whose base rates differ
must have different false positive rates: no score, however it is built,
can equalise all three. The identity is arithmetic, true of any tool on
any data; it says which equalities cannot hold together, not which one
should give way. See `proof("impfpr")`.

## Usage

``` r
impfpr(p, ppv, fnr)
```

## Arguments

- p:

  Base rates, strictly between 0 and 1.

- ppv:

  Positive predictive values, above 0 and at most 1.

- fnr:

  False negative rates, from 0 and below 1. The three have one length
  (or length 1); named ones are matched by name.

## Value

The false positive rates, one per triple.

## References

Chouldechova, A. (2017). Fair prediction with disparate impact: a study
of bias in recidivism prediction instruments. *Big Data*, 5(2), 153-163.
[doi:10.1089/big.2016.0047](https://doi.org/10.1089/big.2016.0047)

## Examples

``` r
# the same ppv and fnr, different base rates: different fpr
impfpr(p = c(a = 0.3, b = 0.5), ppv = 0.6, fnr = 0.25)
#>         a         b 
#> 0.2142857 0.5000000 
```
