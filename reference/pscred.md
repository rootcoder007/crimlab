# Credit for pre-sentence custody

The credit a Canadian court gives for time held in custody before
sentence: each day held counts for at most one day of the sentence
(Criminal Code s. 719(3)) or, where the circumstances justify it, for at
most one and a half days (s. 719(3.1)). With the sentence imposed, the
time left to serve after the credit is given too.

## Usage

``` r
pscred(days, ratio = 1.5, sentence = NULL)
```

## Arguments

- days:

  The days held in custody before sentence: numbers from 0.

- ratio:

  The credit for each day: one number, or one per value of `days`, from
  0 to 1.5.

- sentence:

  `NULL`, or the sentence imposed, in days: one number, or one per value
  of `days`, from 0.

## Value

A data frame with one row per value of `days`: `days`, `ratio`, `credit`
and, with `sentence`, `sentence` and `remain` (the days left to serve).

## Details

The credit is `days * ratio`, in days, unrounded: courts state the
credit they grant (s. 719(3.3)) and round it in their own way. A ratio
of 1.5 is the enhanced credit that R v. Summers (2014) holds is
generally available to an offender denied bail, because time held before
sentence earns no remission and no parole eligibility; a ratio below 1
is a court's choice to give less than the maximum. Any further
mitigation for harsh conditions of custody is a matter of the sentence,
not of this arithmetic. The time left is the sentence less the credit,
and nothing when the credit is the larger.

## References

Criminal Code, R.S.C. 1985, c. C-46, s. 719(3), (3.1) and (3.3).

R. v. Summers, 2014 SCC 26, \[2014\] 1 S.C.R. 575.

## Examples

``` r
pscred(c(10, 30, 45))
#>   days ratio credit
#> 1   10   1.5   15.0
#> 2   30   1.5   45.0
#> 3   45   1.5   67.5
pscred(c(10, 30, 45), ratio = c(1, 1.5, 1.5), sentence = 60)
#>   days ratio credit sentence remain
#> 1   10   1.0   10.0       60     50
#> 2   30   1.5   45.0       60     15
#> 3   45   1.5   67.5       60      0
```
