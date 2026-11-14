# Likelihood ratio of measured evidence

Weighs measurements (a tissue density, a pXRF calcium level, a marker of
antemortem ingestion) under two competing hypotheses, each a normal
model of what the measurement would be: the likelihood ratio is the
probability of the evidence under H1 over that under H2. Its verbal band
follows the usual forensic scale on `|log10 LR|`: 0 none, up to 1 weak,
2 moderate, 3 moderately strong, 4 strong, 6 very strong, above
extremely strong support. The ratio measures support, not the
probability that a hypothesis is true.

## Usage

``` r
preslr(x, m1, s1, m2, s2)
```

## Arguments

- x:

  The measurements: finite numbers.

- m1, s1:

  The mean and standard deviation of a measurement under H1: one number
  each or one per measurement (`s1` positive).

- m2, s2:

  The same under H2.

## Value

A list: `ll1`, `ll2` the log-likelihoods, `llr` their difference (the
natural log of the ratio), `lr` the ratio, `l10` its log10, and
`verbal`, the band and the hypothesis it favours.

## References

ENFSI (2015). *ENFSI Guideline for Evaluative Reporting in Forensic
Science*, approved version 3.0, 8 March 2015. European Network of
Forensic Science Institutes.
<https://enfsi.eu/wp-content/uploads/2016/09/m1_guideline.pdf>

## Examples

``` r
preslr(c(1.9, 2.3, 2.1), m1 = 2, s1 = 0.3, m2 = 1, s2 = 0.4)
#> $ll1
#> [1] 0.2439917
#> 
#> $ll2
#> [1] -11.60169
#> 
#> $llr
#> [1] 11.84569
#> 
#> $lr
#> [1] 139481.2
#> 
#> $l10
#> [1] 5.144516
#> 
#> $verbal
#> [1] "very strong support for H1"
#> 
```
