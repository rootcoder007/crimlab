# Changelog

## crimlab 0.0.0.9000

- [`etasf()`](https://rootcoder007.github.io/crimlab/reference/etasf.md)
  fits the temporal ETAS model (marks and the modified Omori law),
  [`mvhwk()`](https://rootcoder007.github.io/crimlab/reference/mvhwk.md)
  the multivariate Hawkes process of several types that excite each
  other, with its excitation matrix and spectral radius, and
  [`sthwk()`](https://rootcoder007.github.io/crimlab/reference/sthwk.md)
  the space-time Hawkes process of crime, with intensities at query
  points;
  [`hwksim()`](https://rootcoder007.github.io/crimlab/reference/hwksim.md)
  simulates from a space-time fit. `proof("etasf")` names the theorem
  that ties the Omori kernel to the Lomax kernel of
  [`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md).
- [`sdbsim()`](https://rootcoder007.github.io/crimlab/reference/sdbsim.md)
  simulates the Short-D’Orsogna-Brantingham model of burglary hot spots,
  as agents on a lattice or as its reaction-diffusion limit with
  suppressed sites, and
  [`sdbthr()`](https://rootcoder007.github.io/crimlab/reference/sdbthr.md)
  gives the linear stability of its uniform state: the threshold beyond
  which hot spots form and the fastest-growing wavelength.
  `proof("sdbthr")` and `proof("sdbsim")` name the Lean theorems on the
  steady state, the trace and the decay.
- [`hypq()`](https://rootcoder007.github.io/crimlab/reference/hypq.md)
  solves Larson’s hypercube queue of patrol units exactly (workloads,
  dispatch shares, lost or waiting calls, travel times) and
  [`sqrtlw()`](https://rootcoder007.github.io/crimlab/reference/sqrtlw.md)
  gives Kolesar and Blum’s square-root law of response distance;
  `proof("sqrtlw")` names its scaling theorems.
- [`kmrec()`](https://rootcoder007.github.io/crimlab/reference/kmrec.md)
  gives Kaplan-Meier curves by group with Greenwood’s standard errors,
  `survfit()`’s log-scale intervals and medians, and the log-rank test;
  with `unit` it takes event times and studies the gaps between a unit’s
  events (repeat incidents at a place), the last one censored.
- [`trndrr()`](https://rootcoder007.github.io/crimlab/reference/trndrr.md)
  measures the trend in a series of counts as a Poisson rate ratio per
  period (quasi-Poisson when overdispersed), by Mann-Kendall’s test
  (exact up to 8 periods) and the Theil-Sen slope, with the change on
  the period before and a test of equal seasonal totals;
  [`stpcp()`](https://rootcoder007.github.io/crimlab/reference/stpcp.md)
  finds a step change with its permutation p-value and Pettitt’s test;
  [`cusrv()`](https://rootcoder007.github.io/crimlab/reference/cusrv.md)
  watches counts by CUSUM;
  [`xcser()`](https://rootcoder007.github.io/crimlab/reference/xcser.md)
  correlates two series, matched by period name, at lags as
  [`ccf()`](https://rdrr.io/r/stats/acf.html) does.
- [`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md)
  fits a Hawkes (self-exciting) process to event times by maximum
  likelihood with the analytic gradient: exponential, Weibull or Lomax
  kernels, a constant baseline or a trend with a cycle, several of them
  compared by AIC, and the branching arithmetic (cluster size,
  stationary rate, share triggered) that `proof("hwkfit")` names.
  [`hwkres()`](https://rootcoder007.github.io/crimlab/reference/hwkres.md)
  gives the time-rescaling residuals with the Kolmogorov-Smirnov test of
  [`ks.test()`](https://rdrr.io/r/stats/ks.test.html), and
  [`hwksim()`](https://rootcoder007.github.io/crimlab/reference/hwksim.md)
  simulates the process, in time or in space and time, by its branching
  construction.
- [`crsev()`](https://rootcoder007.github.io/crimlab/reference/crsev.md)
  builds a severity-weighted crime index from counts by offence (the
  Crime Severity Index of Statistics Canada: weighted counts per
  population against a base, the weights given or made from
  incarceration rates and mean sentences);
  [`zcomp()`](https://rootcoder007.github.io/crimlab/reference/zcomp.md)
  ranks units on a composite of standardised counts across categories.
- [`cconc()`](https://rootcoder007.github.io/crimlab/reference/cconc.md)
  measures how concentrated a total is across units: the Gini
  coefficient and its exact split into the share of units at zero and
  the Gini of the rest, the Lorenz curve, top shares, dispersion against
  the Poisson, the bound chance alone allows and a simulated Poisson
  null, the Gini within groups, the urban-scaling exponent, and for the
  places of events in order the growth of distinct places.
  [`hilltl()`](https://rootcoder007.github.io/crimlab/reference/hilltl.md)
  estimates a power-law tail index: the discrete maximum likelihood
  estimate through the Hurwitz zeta function, or the Hill estimators,
  with the Kolmogorov-Smirnov distance of the fit. `proof("cconc")`
  names the Lean theorems on the split at zero.
- [`crcap()`](https://rootcoder007.github.io/crimlab/reference/crcap.md)
  estimates a hidden population from two lists (Lincoln-Petersen with
  its bounds over a box of list dependence, Chapman’s estimate with
  Seber’s standard error) or three (the never-listed cell under no
  three-way interaction, its sensitivity to that interaction, and the
  pairwise estimates);
  [`dkfig()`](https://rootcoder007.github.io/crimlab/reference/dkfig.md)
  bounds the dark figure of crime: the true victimisation rate under
  limits on a survey’s false positives and negatives with the breakdown
  point of a conclusion, the total when a survey covers only some
  offences, and the offences behind incidents counted under the
  hierarchy rule. `proof("crcap")` and `proof("dkfig")` name the Lean
  theorems.
- [`bandx()`](https://rootcoder007.github.io/crimlab/reference/bandx.md)
  reads band labels (“1-5”, “Over 60”, “Under 10”) into bounds, places a
  value in each band, expands a banded frequency table and shows how
  much a statistic depends on the cap assumed for an open band;
  [`pbbnd()`](https://rootcoder007.github.io/crimlab/reference/pbbnd.md)
  gives the bounds a rounded or suppressed published cell allows and the
  envelope of the change between two published columns;
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  gives year-over-year (or any period) change within groups with the
  exact interval of the percent change, exact p-values adjusted for
  multiplicity, and the envelope that rounding allows. `proof("pbbnd")`
  names the Lean theorems on the envelope.
- [`jdgiv()`](https://rootcoder007.github.io/crimlab/reference/jdgiv.md)
  gives the judge-leniency instrument’s first stage, intention-to- treat
  effect and Wald ratio with the compliers’ effect over a grid of defier
  shares and effects, or the exact decomposition by type on a known
  population;
  [`jmono()`](https://rootcoder007.github.io/crimlab/reference/jmono.md)
  gives each judge’s decision rate with a chi-square test of equal
  rates, the many-judge slope test of nested judges, and the comparison
  of judges within groups of cases with adjusted p-values and slopes;
  [`mmshr()`](https://rootcoder007.github.io/crimlab/reference/mmshr.md)
  gives the shares of sentences at, below and above a mandatory minimum.
  `proof("jdgiv")` and `proof("jmono")` name the Lean theorems.
- [`agemix()`](https://rootcoder007.github.io/crimlab/reference/agemix.md)
  gives the aggregate age-crime curve of a mixture of latent types and
  the covariance across types between ages that a single type forces to
  0, with its permutation test on panel data;
  [`cumrsk()`](https://rootcoder007.github.io/crimlab/reference/cumrsk.md)
  gives the cumulative risk of a first arrest, conviction or admission
  by age from a life table of first events, with Greenwood’s interval or
  with deaths competing;
  [`jpipe()`](https://rootcoder007.github.io/crimlab/reference/jpipe.md)
  gives rates between stages of the justice pipeline from window sums
  lagged by the time between the stages. `proof("agemix")` names the
  Lean theorem.
- [`mndla()`](https://rootcoder007.github.io/crimlab/reference/mndla.md)
  classifies segregation placements under the Mandela Rules (more than
  15 days, or Rules 43 and 44 with hours out of cell) and gives the
  qualifying share by group over placements, people with any qualifying
  placement and people by their cumulative days, with exact intervals.
- [`stkflw()`](https://rootcoder007.github.io/crimlab/reference/stkflw.md)
  relates person-days, the average daily population, admissions and the
  length of stay (any two give the rest) with changes, rates and the
  stock an earlier release removes;
  [`lostay()`](https://rootcoder007.github.io/crimlab/reference/lostay.md)
  summarises lengths of stay, their bands and thresholds and the
  lognormal, exponential and Pareto fits;
  [`bklog()`](https://rootcoder007.github.io/crimlab/reference/bklog.md)
  applies Little’s law to a docket;
  [`dktflw()`](https://rootcoder007.github.io/crimlab/reference/dktflw.md)
  gives clearance, rolling clearance, the pending caseload, Little’s law
  on published figures and the trend of the outcome mix;
  [`lftab()`](https://rootcoder007.github.io/crimlab/reference/lftab.md)
  builds the actuarial life table;
  [`trnmat()`](https://rootcoder007.github.io/crimlab/reference/trnmat.md)
  gives the transition matrix between states with its stationary
  distribution. `proof("stkflw")` and `proof("bklog")` name the Lean
  theorems.
- [`cntreg()`](https://rootcoder007.github.io/crimlab/reference/cntreg.md)
  fits a count on a group by Poisson and negative binomial regression
  (the iterations of [`glm()`](https://rdrr.io/r/stats/glm.html) and
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)), with
  adjustment columns, an exposure offset and clusters, and reports the
  rate ratios with Wald intervals: model-based, cluster-robust (the
  sandwich of
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)),
  and from an exchangeable GEE (as ‘geepack’ estimates it).
  `proof("cntreg")` names the Lean theorems on swapping the groups.
- [`grprr()`](https://rootcoder007.github.io/crimlab/reference/grprr.md)
  gives rates or proportions by group and their ratios and differences
  against a reference group: Wilson, Clopper-Pearson and exact Poisson
  intervals, log-normal and exact conditional intervals of the ratio,
  Newcombe’s interval of a difference of proportions, and each group’s
  over-representation among the events.
  [`dastd()`](https://rootcoder007.github.io/crimlab/reference/dastd.md)
  gives directly standardised rates with Fay and Feuer’s interval, and
  their ratios and differences;
  [`sirex()`](https://rootcoder007.github.io/crimlab/reference/sirex.md)
  gives standardised incidence ratios with exact intervals, funnel
  limits and Marshall’s empirical Bayes shrinkage.
- [`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
  cross-tabulates labels (counting rows, summing a count column, taking
  the largest or smallest value, or counting distinct ids per cell) and
  tests the table: Pearson’s chi-square (with Yates’ correction on
  request), the likelihood-ratio test, Fisher’s exact test with the
  conditional odds ratio of a 2 x 2 table, Cramer’s V with its
  noncentral chi-square interval and Bergsma’s bias correction, the
  contingency coefficient, phi and the diagonal share. `proof("ctabx")`
  names the Lean theorems on the 2 x 2 statistic.
- [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  maps a published dataset’s columns onto a generic function’s arguments
  from `inst/presets/presets.tsv` (one row per preset and argument:
  column, recode, filter, units, note, source), with recipes for derived
  columns, so that `generic(clprep(data, preset))` gives the published
  analysis. The presets cover the Ontario segregation and
  deaths-in-custody tables, the race-based use-of-force files, the
  Special Investigations Unit cases, the Structured Intervention Unit
  tables and the Toronto incident layers. Presets also feed the
  risk-tool audits:
  [`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md),
  [`tmpaud()`](https://rootcoder007.github.io/crimlab/reference/tmpaud.md),
  [`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md)
  and
  [`scrdsp()`](https://rootcoder007.github.io/crimlab/reference/scrdsp.md)
  take a
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  result as their first argument, for a generic file of area, risk,
  outcome and group, for
  [`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md)
  records, and for the layouts of two published studies of risk
  classification.
- Audits of risk tools and predictive policing.
  [`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md)
  gives the parity gap and disparate impact ratio between groups with
  Newcombe and Katz intervals, the gaps in true and false positive rates
  and the average odds difference with their intervals, and the bias
  amplification score with a bootstrap interval;
  [`tmpaud()`](https://rootcoder007.github.io/crimlab/reference/tmpaud.md)
  the same by city and period with each city’s range.
  [`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md)
  rolls records up to areas and asks whether the order of predicted risk
  matches the order of recorded outcomes (Spearman’s rho with the AS 89
  and large-sample p-values) and whether the mean rank gap of each
  group’s areas departs from chance (a permutation test);
  [`scrdsp()`](https://rootcoder007.github.io/crimlab/reference/scrdsp.md)
  a score’s spread across groups with the one-way analysis of variance
  and Welch’s test.
  [`lipaud()`](https://rootcoder007.github.io/crimlab/reference/lipaud.md)
  is the Lipschitz audit of individual fairness.
  [`nsbnd()`](https://rootcoder007.github.io/crimlab/reference/nsbnd.md)
  bounds a true base rate behind a noisily recorded one for noise rates
  in boxes, and says which groups’ order survives every admissible
  noise;
  [`impfpr()`](https://rootcoder007.github.io/crimlab/reference/impfpr.md)
  is Chouldechova’s identity;
  [`rnkres()`](https://rootcoder007.github.io/crimlab/reference/rnkres.md)
  the pairs a noisy ranking orders.
  [`oaxbd()`](https://rootcoder007.github.io/crimlab/reference/oaxbd.md)
  is the Oaxaca-Blinder decomposition of a gap at either group’s
  coefficients and the threefold one, with delta-method standard errors.
  [`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md)
  simulates records with a known disparity and
  [`nordet()`](https://rootcoder007.github.io/crimlab/reference/nordet.md)
  the noisy-OR chance of detection by patrols.
  [`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
  names the Lean theorems behind
  [`nsbnd()`](https://rootcoder007.github.io/crimlab/reference/nsbnd.md),
  [`impfpr()`](https://rootcoder007.github.io/crimlab/reference/impfpr.md),
  [`rnkres()`](https://rootcoder007.github.io/crimlab/reference/rnkres.md)
  and
  [`oaxbd()`](https://rootcoder007.github.io/crimlab/reference/oaxbd.md).
- [`relchk()`](https://rootcoder007.github.io/crimlab/reference/relchk.md)
  recomputes a reported odds ratio, risk ratio or risk difference from
  its table under every labelling (each reference group, each outcome as
  the event, groups in rows or in columns) and says which readings
  reproduce the value as rounded: the labelling the report states,
  another one, or none. Zero cells give `Inf` or `NaN`, never a
  correction. `proof("relchk")` names the Lean theorems on how the
  contrasts change under relabelling.
- [`catchk()`](https://rootcoder007.github.io/crimlab/reference/catchk.md)
  audits categorical inputs before an analysis: storage type, levels and
  counts, missing values, unused and rare levels, labels that differ
  only by case or white space, 0/1 and logical codings, numbers that
  look like category codes, value labels that disagree with the values,
  and the reference level a default would take and why. It can check a
  stated reference level, and changes nothing.
- Disparity tests:
  [`vodtst()`](https://rootcoder007.github.io/crimlab/reference/vodtst.md),
  [`orrisk()`](https://rootcoder007.github.io/crimlab/reference/orrisk.md),
  [`outtst()`](https://rootcoder007.github.io/crimlab/reference/outtst.md),
  [`dspbm()`](https://rootcoder007.github.io/crimlab/reference/dspbm.md),
  [`dyadrt()`](https://rootcoder007.github.io/crimlab/reference/dyadrt.md).
- Forensic statistics:
  [`gtmp()`](https://rootcoder007.github.io/crimlab/reference/gtmp.md),
  [`lrpost()`](https://rootcoder007.github.io/crimlab/reference/lrpost.md),
  [`cllr()`](https://rootcoder007.github.io/crimlab/reference/cllr.md),
  [`cllrmin()`](https://rootcoder007.github.io/crimlab/reference/cllrmin.md),
  [`lrcal()`](https://rootcoder007.github.io/crimlab/reference/lrcal.md),
  [`mixlr()`](https://rootcoder007.github.io/crimlab/reference/mixlr.md),
  [`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md),
  [`patidx()`](https://rootcoder007.github.io/crimlab/reference/patidx.md),
  [`dbsrch()`](https://rootcoder007.github.io/crimlab/reference/dbsrch.md),
  [`tipplt()`](https://rootcoder007.github.io/crimlab/reference/tipplt.md),
  [`eceplt()`](https://rootcoder007.github.io/crimlab/reference/eceplt.md).
- [`clrsmp()`](https://rootcoder007.github.io/crimlab/reference/clrsmp.md)
  splits a change in a clearance rate into rate and mix effects and
  flags Simpson’s paradox.
- [`expbnd()`](https://rootcoder007.github.io/crimlab/reference/expbnd.md)
  bounds rates and disparities when exposure is known only through a
  proxy.
- [`misrd()`](https://rootcoder007.github.io/crimlab/reference/misrd.md)
  gives the risk difference when group membership is misrecorded.
- [`rtmhs()`](https://rootcoder007.github.io/crimlab/reference/rtmhs.md)
  measures regression to the mean at selected hot spots.
- [`incap()`](https://rootcoder007.github.io/crimlab/reference/incap.md)
  gives the crime rate and share prevented under incapacitation.
- [`incapr()`](https://rootcoder007.github.io/crimlab/reference/incapr.md)
  gives crimes prevented under desistance and replacement.
- [`sellab()`](https://rootcoder007.github.io/crimlab/reference/sellab.md)
  bounds the failure rate of a release rule from released defendants
  alone.
- [`losbnd()`](https://rootcoder007.github.io/crimlab/reference/losbnd.md)
  bounds the mean time in custody from completed and open stays.
- [`sntbnd()`](https://rootcoder007.github.io/crimlab/reference/sntbnd.md)
  bounds the effect of a sentence on a binary outcome without
  assumptions.
- [`hzdsel()`](https://rootcoder007.github.io/crimlab/reference/hzdsel.md)
  shows how period-two hazards differ after selective depletion.
- Networks:
  [`coofg()`](https://rootcoder007.github.io/crimlab/reference/coofg.md)
  builds the co-offending graph of an incident table;
  [`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md),
  [`ntbtw()`](https://rootcoder007.github.io/crimlab/reference/ntbtw.md),
  [`ntcls()`](https://rootcoder007.github.io/crimlab/reference/ntcls.md)
  and
  [`ntevc()`](https://rootcoder007.github.io/crimlab/reference/ntevc.md)
  give degree, betweenness, closeness and eigenvector centrality;
  [`ntcore()`](https://rootcoder007.github.io/crimlab/reference/ntcore.md)
  cores,
  [`ntcmp()`](https://rootcoder007.github.io/crimlab/reference/ntcmp.md)
  components,
  [`ntlvn()`](https://rootcoder007.github.io/crimlab/reference/ntlvn.md)
  and
  [`ntmod()`](https://rootcoder007.github.io/crimlab/reference/ntmod.md)
  Louvain communities and modularity,
  [`ntbrk()`](https://rootcoder007.github.io/crimlab/reference/ntbrk.md)
  Gould-Fernandez brokerage and
  [`ntchg()`](https://rootcoder007.github.io/crimlab/reference/ntchg.md)
  the Cheeger bounds. Results are by node label in order of first
  appearance.
- Spatial:
  [`hvdist()`](https://rootcoder007.github.io/crimlab/reference/hvdist.md)
  and
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)/[`eqinv()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  for longitude and latitude;
  [`cgstat()`](https://rootcoder007.github.io/crimlab/reference/cgstat.md),
  [`cgmedn()`](https://rootcoder007.github.io/crimlab/reference/cgmedn.md),
  [`webpt()`](https://rootcoder007.github.io/crimlab/reference/webpt.md),
  [`sdell()`](https://rootcoder007.github.io/crimlab/reference/sdell.md)
  and
  [`cvhull()`](https://rootcoder007.github.io/crimlab/reference/cvhull.md)
  for centres, spread and hulls;
  [`nnidx()`](https://rootcoder007.github.io/crimlab/reference/nnidx.md),
  [`ripk()`](https://rootcoder007.github.io/crimlab/reference/ripk.md)
  and
  [`ripenv()`](https://rootcoder007.github.io/crimlab/reference/ripenv.md)
  for nearest neighbours and K functions;
  [`fzmode()`](https://rootcoder007.github.io/crimlab/reference/fzmode.md),
  [`kmpart()`](https://rootcoder007.github.io/crimlab/reference/kmpart.md),
  [`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)
  and
  [`stacel()`](https://rootcoder007.github.io/crimlab/reference/stacel.md)
  for hot spots;
  [`mantst()`](https://rootcoder007.github.io/crimlab/reference/mantst.md)
  and
  [`stmavg()`](https://rootcoder007.github.io/crimlab/reference/stmavg.md)
  in space and time;
  [`jtdecay()`](https://rootcoder007.github.io/crimlab/reference/jtdecay.md),
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md),
  [`jtcal()`](https://rootcoder007.github.io/crimlab/reference/jtcal.md),
  [`jtacc()`](https://rootcoder007.github.io/crimlab/reference/jtacc.md)
  and
  [`bjtc()`](https://rootcoder007.github.io/crimlab/reference/bjtc.md)
  for journey to crime. Points come as vectors, tables or sf points;
  longitude and latitude only with `lonlat = TRUE`.
- Spatial weights and autocorrelation:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)
  (k nearest, distance band, contiguity; binary or inverse distance),
  [`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md)
  (Geary’s C),
  [`gorg()`](https://rootcoder007.github.io/crimlab/reference/gorg.md)
  (the Getis-Ord G),
  [`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md)
  (local Moran with conditional permutations),
  [`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md)
  (local G and G\*) and
  [`dirmn()`](https://rootcoder007.github.io/crimlab/reference/dirmn.md)
  (directional mean).
- [`hdbang()`](https://rootcoder007.github.io/crimlab/reference/hdbang.md)
  (head-bang smoothing, weighted or not),
  [`crwalk()`](https://rootcoder007.github.io/crimlab/reference/crwalk.md)
  (a correlated random walk of a series: fit, expected displacement and
  simulation) and
  [`rnnhcl()`](https://rootcoder007.github.io/crimlab/reference/rnnhcl.md)
  (nearest-neighbour clustering against a population at risk).
  [`nnidx()`](https://rootcoder007.github.io/crimlab/reference/nnidx.md),
  [`ripk()`](https://rootcoder007.github.io/crimlab/reference/ripk.md),
  [`fzmode()`](https://rootcoder007.github.io/crimlab/reference/fzmode.md)
  and
  [`mantst()`](https://rootcoder007.github.io/crimlab/reference/mantst.md)
  take `metric = "manhattan"`;
  [`fzmode()`](https://rootcoder007.github.io/crimlab/reference/fzmode.md),
  [`kmpart()`](https://rootcoder007.github.io/crimlab/reference/kmpart.md),
  [`nnhclu()`](https://rootcoder007.github.io/crimlab/reference/nnhclu.md)
  and
  [`stacel()`](https://rootcoder007.github.io/crimlab/reference/stacel.md)
  take zones as counts at centroids (`count`).
- [`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md)
  (fairness of a risk tool by group),
  [`tlvc()`](https://rootcoder007.github.io/crimlab/reference/tlvc.md)
  and
  [`tlvlr()`](https://rootcoder007.github.io/crimlab/reference/tlvc.md)
  (trace evidence),
  [`aorwt()`](https://rootcoder007.github.io/crimlab/reference/aorwt.md)
  (aoristic weights),
  [`bdig()`](https://rootcoder007.github.io/crimlab/reference/bdig.md)
  and
  [`bnft()`](https://rootcoder007.github.io/crimlab/reference/bdig.md)
  (Benford screens),
  [`skp()`](https://rootcoder007.github.io/crimlab/reference/skp.md)
  (Scott-Knott grouping),
  [`rptc()`](https://rootcoder007.github.io/crimlab/reference/rptc.md)
  (constant reporters),
  [`rswp()`](https://rootcoder007.github.io/crimlab/reference/rswp.md)
  and
  [`rshl()`](https://rootcoder007.github.io/crimlab/reference/rswp.md)
  (Rosenbaum bounds),
  [`zgll()`](https://rootcoder007.github.io/crimlab/reference/zgll.md),
  [`zgem()`](https://rootcoder007.github.io/crimlab/reference/zgll.md),
  [`mvtll()`](https://rootcoder007.github.io/crimlab/reference/mvtll.md)
  and
  [`mvtem()`](https://rootcoder007.github.io/crimlab/reference/mvtll.md)
  (group-based trajectories),
  [`tfrp()`](https://rootcoder007.github.io/crimlab/reference/tfrp.md)
  (glass transfer),
  [`lwss()`](https://rootcoder007.github.io/crimlab/reference/lwss.md)
  (lowess),
  [`chmb()`](https://rootcoder007.github.io/crimlab/reference/chmb.md)
  (tool marks),
  [`gstr()`](https://rootcoder007.github.io/crimlab/reference/gstr.md)
  (grid Gi\*) and
  [`sptst()`](https://rootcoder007.github.io/crimlab/reference/sptst.md)
  (Andresen’s point pattern test).
- [`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
  names the Lean theorem behind each function; the proofs ship in
  `inst/lean`.
- A vector longer than the largest integer is refused with an error, by
  the exported functions and by the compiled entry points alike; the
  compiled code used to truncate its length.
- [`eceplt()`](https://rootcoder007.github.io/crimlab/reference/eceplt.md)
  gives the calibrated curve its limits at extreme prior odds: a
  calibrated ratio of 0 or infinity used to give `NaN`, or an error.
- [`incapr()`](https://rootcoder007.github.io/crimlab/reference/incapr.md)
  refuses a `t0 + s` beyond the length of `lam` however large the two
  are; past the largest integer it used to warn of an integer overflow
  and stop with an unrelated error.
- [`sntbnd()`](https://rootcoder007.github.io/crimlab/reference/sntbnd.md)
  orders the two values of `z` by Unicode code point, so the default
  `treat` and the order of the terms are the same on every machine and
  in every locale (`"B"` sorts before `"a"`).
- Every error message starts with the name of the function that raised
  it (`` "outtst: `q` must be at most 1." ``); the words after it are
  unchanged.
- A matrix with more than one row and column is refused where a vector
  is expected, naming the argument; it used to be read as one long
  vector.
- [`eceplt()`](https://rootcoder007.github.io/crimlab/reference/eceplt.md)
  labels its vertical axis in bits.
- [`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
  names the theorems now proved for
  [`gtmp()`](https://rootcoder007.github.io/crimlab/reference/gtmp.md),
  [`cllr()`](https://rootcoder007.github.io/crimlab/reference/cllr.md),
  [`dbsrch()`](https://rootcoder007.github.io/crimlab/reference/dbsrch.md),
  [`mixlr()`](https://rootcoder007.github.io/crimlab/reference/mixlr.md),
  [`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md),
  [`patidx()`](https://rootcoder007.github.io/crimlab/reference/patidx.md),
  [`misrd()`](https://rootcoder007.github.io/crimlab/reference/misrd.md),
  [`incap()`](https://rootcoder007.github.io/crimlab/reference/incap.md),
  [`lrcal()`](https://rootcoder007.github.io/crimlab/reference/lrcal.md)
  and
  [`cllrmin()`](https://rootcoder007.github.io/crimlab/reference/cllrmin.md).
  The help pages of
  [`misrd()`](https://rootcoder007.github.io/crimlab/reference/misrd.md),
  [`cllr()`](https://rootcoder007.github.io/crimlab/reference/cllr.md)
  and
  [`cllrmin()`](https://rootcoder007.github.io/crimlab/reference/cllrmin.md)
  say what is proved and what is not.
- [`vignette("crimlab")`](https://rootcoder007.github.io/crimlab/articles/crimlab.md)
  walks through one example per family.
- Arguments that hold one value per group are matched by name. Named
  inputs whose names come in another order were read by position, a
  wrong answer:
  `vodtst(c(white = 600, black = 400), c(black = 300, white = 700))`
  read 300 as the white night count. Now every later argument is put in
  the order of the first by name; names that differ, or a name given
  twice, are an error that lists them. Unnamed inputs are still matched
  by position. This holds in
  [`vodtst()`](https://rootcoder007.github.io/crimlab/reference/vodtst.md),
  [`dspbm()`](https://rootcoder007.github.io/crimlab/reference/dspbm.md),
  [`dyadrt()`](https://rootcoder007.github.io/crimlab/reference/dyadrt.md)
  (the row and column names of `count` too),
  [`clrsmp()`](https://rootcoder007.github.io/crimlab/reference/clrsmp.md),
  [`expbnd()`](https://rootcoder007.github.io/crimlab/reference/expbnd.md),
  [`incap()`](https://rootcoder007.github.io/crimlab/reference/incap.md),
  [`gtmp()`](https://rootcoder007.github.io/crimlab/reference/gtmp.md),
  [`hzdsel()`](https://rootcoder007.github.io/crimlab/reference/hzdsel.md),
  and for loci in
  [`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md)
  and
  [`patidx()`](https://rootcoder007.github.io/crimlab/reference/patidx.md);
  per-unit arguments of
  [`outtst()`](https://rootcoder007.github.io/crimlab/reference/outtst.md),
  [`rtmhs()`](https://rootcoder007.github.io/crimlab/reference/rtmhs.md),
  [`lrpost()`](https://rootcoder007.github.io/crimlab/reference/lrpost.md),
  [`sellab()`](https://rootcoder007.github.io/crimlab/reference/sellab.md)
  and
  [`sntbnd()`](https://rootcoder007.github.io/crimlab/reference/sntbnd.md)
  are matched by name when every one has unique names.
- Per-group arguments also accept a list of single numbers and a one-row
  data frame; a one-row or one-column matrix gives its names.
- [`vodtst()`](https://rootcoder007.github.io/crimlab/reference/vodtst.md)
  takes one table of groups by day and night (a matrix, 2-d table or
  data frame with columns named `day` and `night`) with `night` left
  out.
- [`vodtst()`](https://rootcoder007.github.io/crimlab/reference/vodtst.md),
  [`dspbm()`](https://rootcoder007.github.io/crimlab/reference/dspbm.md)
  and
  [`expbnd()`](https://rootcoder007.github.io/crimlab/reference/expbnd.md)
  return a `ref` column marking the reference group; a factor `ref` is
  read by its label.
- [`sntbnd()`](https://rootcoder007.github.io/crimlab/reference/sntbnd.md)
  no longer picks the treatment by sorting the values of `z`: `treat` is
  required unless `z` is logical (`TRUE` is the treatment) or 0 and 1 (1
  is).
- A labelled vector (with a `labels` attribute, as haven makes) is
  refused with a request to decode it first: its numbers are category
  codes.
- Learning:
  [`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md),
  the elastic net path by coordinate descent with cross-validation over
  given folds (the objective and scale of glmnet);
  [`lgtreg()`](https://rootcoder007.github.io/crimlab/reference/lgtreg.md),
  logistic regression;
  [`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md),
  random forests for regression and probabilities.
- Causal estimates:
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
  double machine learning by cross-fitting for the partially linear,
  interactive, IV, difference-in-differences and sample-selection
  models, with AIPW and TMLE from the same nuisances.
- [`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md)
  screens binomial paths by the sequential strong rule (`scr`);
  [`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md)
  grows honest trees (`hon`), limits split imbalance (`imb`) and takes
  `mtry` by the rule of grf by default;
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
  adds two-way clusters, the experimental DiD scores (`score = "exp"`)
  and passes `tol`, `scr` and `hon` to its learners.
- Policy evaluation on unit-by-period panels.
  [`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
  is the two-way fixed-effects difference in differences by least
  squares (alternating projections, weights, one-way or two-way
  cluster-robust errors with the small-sample factor of ‘fixest’), or by
  Poisson or negative binomial regression for counts, with the estimate
  leaving out each treated unit, the Goodman-Bacon decomposition, and
  matching on the outcomes before treatment.
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  is the event study with the joint test of the leads and the raw trends
  of each adoption cohort.
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md)
  gives the group-time effects of Callaway and Sant’Anna, with the
  doubly robust, IPW and regression estimators for covariates,
  Wooldridge’s regression on cohort-period cells, and de Chaisemartin
  and D’Haultfoeuille’s estimator for treatments that switch on and off,
  with cluster-bootstrap errors.
  [`synthc()`](https://rootcoder007.github.io/crimlab/reference/synthc.md)
  is the synthetic control of one treated unit with the placebo test,
  and
  [`sdidx()`](https://rootcoder007.github.io/crimlab/reference/sdidx.md)
  the synthetic difference in differences with placebo or jackknife
  errors; both solve their weights exactly (Wolfe’s minimum-norm point).
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  bounds the effects under departures from parallel trends (Rambachan
  and Roth’s relative magnitudes) with a breakdown value.
  [`wcbx()`](https://rootcoder007.github.io/crimlab/reference/wcbx.md)
  gives the restricted wild cluster bootstrap (every sign pattern when
  few), permutation of adoption dates across clusters, and simulated
  power with the minimum detectable effect;
  [`wystep()`](https://rootcoder007.github.io/crimlab/reference/wystep.md)
  the Westfall-Young stepdown p-values.
- Spatial crime analysis.
  [`knxnr()`](https://rootcoder007.github.io/crimlab/reference/knxnr.md)
  is the Knox test of near repeats: pairs by band of distance and time
  against their exact expectation under exchangeable times, with
  permutation and Poisson p-values.
  [`kdsurf()`](https://rootcoder007.github.io/crimlab/reference/kdsurf.md)
  gives kernel intensity surfaces (the bandwidth and values of
  [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html), the edge
  corrections of
  [`spatstat.explore::density.ppp()`](https://rdrr.io/pkg/spatstat.explore/man/density.ppp.html))
  with the smallest areas holding each share of the incidents;
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md)
  DBSCAN clusters as `dbscan::dbscan()` numbers them;
  [`hsper()`](https://rootcoder007.github.io/crimlab/reference/hsper.md)
  hot grid cells per period, their overlap and how long each stays hot;
  [`ebshr()`](https://rootcoder007.github.io/crimlab/reference/ebshr.md)
  the empirical-Bayes shrinkage of hot-spot counts and the fall to
  expect without any programme;
  [`kscan()`](https://rootcoder007.github.io/crimlab/reference/kscan.md)
  Kulldorff’s space-time permutation scan with secondary clusters.
  [`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)
  and
  [`bvmor()`](https://rootcoder007.github.io/crimlab/reference/bvmor.md)
  give global and bivariate Moran’s I with their inference
  ([`spdep::moran.test()`](https://r-spatial.github.io/spdep/reference/moran.test.html),
  [`spdep::moran_bv()`](https://r-spatial.github.io/spdep/reference/moran_bv.html)),
  one variable or a sweep of columns;
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  fits the spatial lag, error and combined models by maximum likelihood,
  and the lag model by two-stage least squares, with Durbin terms and
  the OLS fit beside (`spatialreg::lagsarlm()`, `errorsarlm()`,
  `sacsarlm()`, `stsls()`);
  [`spill()`](https://rootcoder007.github.io/crimlab/reference/spill.md)
  gives direct, indirect and total effects with delta-method standard
  errors (`spatialreg::impacts()`).
  [`crossk()`](https://rootcoder007.github.io/crimlab/reference/crossk.md)
  gives the cross-type K function and the kernel segregation test,
  [`qdrtst()`](https://rootcoder007.github.io/crimlab/reference/qdrtst.md)
  the quadrat test,
  [`morsta()`](https://rootcoder007.github.io/crimlab/reference/morsta.md)
  Morisita’s index,
  [`dissix()`](https://rootcoder007.github.io/crimlab/reference/dissix.md)
  the dissimilarity, entropy, Atkinson and Gini indices with Morrill’s
  adjustment.
  [`pinpl()`](https://rootcoder007.github.io/crimlab/reference/pinpl.md)
  puts points in polygons and reports the population their regions
  cover,
  [`arwt()`](https://rootcoder007.github.io/crimlab/reference/arwt.md)
  moves counts between two sets of areal units, and
  [`bufas()`](https://rootcoder007.github.io/crimlab/reference/bufas.md)
  assigns treated, comparison and excluded places around sites.
  `proof("ebshr")`, `proof("kscan")`, `proof("dissix")`,
  `proof("spill")` and `proof("bvmor")` name the Lean theorems in
  `inst/lean/Hotspot.lean`.
- Design-based causal effects of a 0/1 treatment.
  [`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md)
  gives Hajek-normalised inverse probability weighting, g-computation,
  subclassification and caliper matching on the propensity (ATE, ATT or
  ATC), and entropy balancing for the ATT, with bootstrap errors;
  [`smdbal()`](https://rootcoder007.github.io/crimlab/reference/smdbal.md)
  the standardised mean differences and variance ratios before and after
  weighting with the overlap of the propensity;
  [`evlx()`](https://rootcoder007.github.io/crimlab/reference/evlx.md)
  the E-value of a risk, odds or hazard ratio or a standardised
  difference.
- [`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
  gives design-based means, totals, proportions, ratios and variances
  from stratified, clustered and weighted samples, by domain and with
  differences between domains, with linearisation errors, the
  finite-population correction and a description of the design.
