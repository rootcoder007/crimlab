# Concentration-response functions for PM2.5 and NO2

A concentration-response function (CRF) turns an air-pollution
concentration into the relative risk (RR) of a health outcome compared
with a counterfactual “clean air” level: an RR of 1.05 means 5% more
deaths (or cases) than at the reference concentration. These are the
first step of any burden estimate;
[`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
calls them.

## Usage

``` r
crfpm25(expo, y = "allcause", refconc = 5.8)

crfno2(expo, y = "allcause", refconc = 10, betaper10 = NULL)

# S3 method for class 'crf'
print(x, ...)
```

## Arguments

- expo:

  Ambient concentration in micrograms per cubic metre, a single
  non-negative number or a vector of them (one per area, day or person).

- y:

  PM2.5: `"allcause"`, `"ihd"` or `"stroke"` (deaths). NO2: `"allcause"`
  and `"respiratory"` (deaths; Huangfu and Atkinson 2020, RR 1.02 and
  1.03 per 10 micrograms per cubic metre) or `"chiast"` (incident
  childhood asthma, not deaths; Khreis et al. 2017, OR 1.05 per 4
  micrograms per cubic metre, scaled to per 10).

- refconc:

  Counterfactual concentration below which no excess risk is assumed.
  PM2.5 default 5.8, the lower bound of the GBD 2019 theoretical minimum
  risk exposure level (the WHO 2021 guideline is 5); NO2 default 10, the
  WHO 2021 guideline.

- betaper10:

  Optional log relative risk per 10 micrograms per cubic metre (a single
  finite number, e.g. `log(1.04)`) overriding the NO2 outcome table;
  `NULL` (the default) uses the table.

- x:

  An object of class `crf`.

- ...:

  Not used.

## Value

A list of class `crf` with: `rr` and `logrr`, the relative risk and its
logarithm at the exposure (their means over a vector of exposures);
`rrlo` and `rrup`, the RR at the published 95% limits (`NA` where there
are none); `refconc`; `expconc` (the exposure, or the mean of a vector);
`pollutant` (`"PM2.5"` or `"NO2"`); `citation` (the source as a string);
and `extra`, a list with the coefficient (`betaper10` and `betaper10ci`
for the log-linear forms, `alpha`, `gamma`, `delta` for the IER), `form`
(PM2.5 only), `outcome`, and `rrperunit`, the RR at each exposure
(`NULL` for one exposure).

## Details

They belong to environmental criminology in both its established senses:
the study of environmental crime and harm (green criminology), where the
burden a polluting activity places on an area or group, and how unevenly
it falls, is the measure of the harm; and the study of the environment
as a driver of offending and victimisation, where exposure series are
joined to crime counts by place and day, as in work relating air
pollution to crime (Bondy, Roth and Sager 2020; Herrnstadt et al. 2019).
The second example below is of that kind.

`crfpm25()` is log-linear for all-cause mortality with the pooled cohort
estimate of the WHO 2021 guideline review (Chen and Hoek 2020: RR 1.08,
95% CI 1.06-1.09, per 10 micrograms per cubic metre), \$\$RR(z) =
\exp(\ln(1.08) (z - z\_{cf}) / 10)\$\$ for \\z \> z\_{cf}\\ and 1
otherwise; for the cause-specific outcomes (ischaemic heart disease,
stroke) it is the Integrated Exposure-Response curve of Burnett et al.
(2014, equation 1), \$\$RR(z) = 1 + \alpha (1 - e^{-\gamma (z -
z\_{cf})^{\delta}}),\$\$ with the GBD 2013 parameter triples; the IER
was fitted per cause and has no all-cause form. `crfno2()` is the
log-linear function of the WHO 2021 review (Huangfu and Atkinson 2020:
RR 1.02, 95% CI 1.01-1.04, per 10 micrograms per cubic metre for
all-cause mortality), \$\$RR(z) = \exp(\beta (z - z\_{cf}) / 10).\$\$ A
vector of exposures gives the mean RR and log-RR; the per-unit relative
risks are then in `extra$rrperunit`.

Every function has zero excess risk at and below the reference
concentration (\\z \le z\_{cf}\\ gives RR 1). The published 95% limits
of the log-linear coefficients are carried through the same curve into
`rrlo` and `rrup`; the IER parameters come without an interval, so those
are `NA` for `"ihd"` and `"stroke"`, and also when `betaper10` is
supplied.

Refusals: an unknown `y` (the message lists the available ones); an
exposure that is missing, infinite, negative, logical or empty; and a
`refconc` (or `betaper10`) that is not a single finite number, or a
negative `refconc`. A supplied `betaper10` may be negative (a protective
association) and replaces the table, in which case `y` is only a label.

## References

Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
contemporaneous relationship between air pollution and crime. Journal of
the Association of Environmental and Resource Economists 7(3), 555-585.
[doi:10.1086/707127](https://doi.org/10.1086/707127) . Herrnstadt, E.,
Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air pollution as a
cause of violent crime. Working paper. Chen, J. and Hoek, G. (2020).
Long-term exposure to PM and all-cause and cause-specific mortality: a
systematic review and meta-analysis. Environment International 143,
105974. Huangfu, P. and Atkinson, R. (2020). Long-term exposure to NO2
and O3 and all-cause and respiratory mortality: a systematic review and
meta-analysis. Environment International 144, 105998. Khreis, H. et al.
(2017). Exposure to traffic-related air pollution and risk of
development of childhood asthma: a systematic review and meta-analysis.
Environment International 100, 1-31. GBD 2019 Risk Factors Collaborators
(2020). The Lancet 396(10258), 1223-1249. Burnett, R. T. et al. (2014).
An integrated risk function for estimating the global burden of disease
attributable to ambient fine particulate matter exposure. Environmental
Health Perspectives 122(4), 397-403. WHO (2021). Global Air Quality
Guidelines.

## Examples

``` r
# All-cause mortality at an annual PM2.5 of 12 micrograms per cubic metre
r <- crfpm25(12)
c(rr = r$rr, lower = r$rrlo, upper = r$rrup)
#>       rr    lower    upper 
#> 1.048873 1.036787 1.054883 
# the cause-specific IER curve for ischaemic heart disease
crfpm25(12, y = "ihd")$rr
#> [1] 1.553634
# NO2 and respiratory deaths; at or below the reference no excess risk
crfno2(25, y = "respiratory")$rr
#> [1] 1.045336
crfno2(8)$rr
#> [1] 1
# a vector of exposures (one per area): the mean RR and each area's RR
crfpm25(c(5, 12, 20))$extra$rrperunit
#> [1] 1.000000 1.048873 1.115480
# a coefficient of the reader's own: RR 1.04 per 10 micrograms
crfno2(30, betaper10 = log(1.04))$rr
#> [1] 1.0816
# the environment as a driver of offending: a simulated daily NO2 series
# and simulated daily counts of assaults for three areas, joined by area
# and day
set.seed(1)
d <- data.frame(area = rep(c("a", "b", "c"), each = 30), day = rep(1:30, 3))
d$no2 <- stats::rlnorm(90, log(20), 0.3)
d$assaults <- stats::rpois(90, 3 + 0.05 * d$no2)
# the relative risk each area-day's exposure carries, and the mean count
# on days above and below its median
d$rr <- crfno2(d$no2)$extra$rrperunit
tapply(d$assaults, d$rr > stats::median(d$rr), mean)
#>    FALSE     TRUE 
#> 3.822222 4.222222 
# an outcome without a published function is refused
try(crfpm25(12, y = "asthma"))
#> Error : crfpm25: Unknown outcome 'asthma'. Available: allcause, ihd, stroke.
```
