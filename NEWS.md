# crimlab 0.0.0.9000

* `trndrr()` measures the trend in a series of counts as a Poisson rate
  ratio per period (quasi-Poisson when overdispersed), by Mann-Kendall's
  test (exact up to 8 periods) and the Theil-Sen slope, with the change on
  the period before and a test of equal seasonal totals; `stpcp()` finds a
  step change with its permutation p-value and Pettitt's test; `cusrv()`
  watches counts by CUSUM; `xcser()` correlates two series, matched by
  period name, at lags as `ccf()` does.
* `hwkfit()` fits a Hawkes (self-exciting) process to event times by
  maximum likelihood with the analytic gradient: exponential, Weibull or
  Lomax kernels, a constant baseline or a trend with a cycle, several of
  them compared by AIC, and the branching arithmetic (cluster size,
  stationary rate, share triggered) that `proof("hwkfit")` names.
  `hwkres()` gives the time-rescaling residuals with the Kolmogorov-Smirnov
  test of `ks.test()`, and `hwksim()` simulates the process, in time or in
  space and time, by its branching construction.
* `cntreg()` fits a count on a group by Poisson and negative binomial
  regression (the iterations of `glm()` and `MASS::glm.nb()`), with
  adjustment columns, an exposure offset and clusters, and reports the
  rate ratios with Wald intervals: model-based, cluster-robust (the
  sandwich of `sandwich::vcovCL()`), and from an exchangeable GEE (as
  'geepack' estimates it). `proof("cntreg")` names the Lean theorems on
  swapping the groups.
* `grprr()` gives rates or proportions by group and their ratios and
  differences against a reference group: Wilson, Clopper-Pearson and
  exact Poisson intervals, log-normal and exact conditional intervals of
  the ratio, Newcombe's interval of a difference of proportions, and each
  group's over-representation among the events. `dastd()` gives directly
  standardised rates with Fay and Feuer's interval, and their ratios and
  differences; `sirex()` gives standardised incidence ratios with exact
  intervals, funnel limits and Marshall's empirical Bayes shrinkage.
* `ctabx()` cross-tabulates labels (counting rows, summing a count column,
  taking the largest or smallest value, or counting distinct ids per
  cell) and tests the table: Pearson's chi-square (with Yates' correction
  on request), the likelihood-ratio test, Fisher's exact test with the
  conditional odds ratio of a 2 x 2 table, Cramer's V with its noncentral
  chi-square interval and Bergsma's bias correction, the contingency
  coefficient, phi and the diagonal share. `proof("ctabx")` names the
  Lean theorems on the 2 x 2 statistic.
* `clprep()` maps a published dataset's columns onto a generic function's
  arguments from `inst/presets/presets.tsv` (one row per preset and
  argument: column, recode, filter, units, note, source), with recipes
  for derived columns, so that `generic(clprep(data, preset))` gives the
  published analysis. The presets cover the Ontario segregation and
  deaths-in-custody tables, the race-based use-of-force files, the
  Special Investigations Unit cases, the Structured Intervention Unit
  tables and the Toronto incident layers.
  Presets also feed the risk-tool audits: `fairdf()`, `tmpaud()`,
  `calaud()` and `scrdsp()` take a `clprep()` result as their first
  argument, for a generic file of area, risk, outcome and group, for
  `fairsim()` records, and for the layouts of two published studies of
  risk classification.
* Audits of risk tools and predictive policing. `fairdf()` gives the
  parity gap and disparate impact ratio between groups with Newcombe and
  Katz intervals, the gaps in true and false positive rates and the
  average odds difference with their intervals, and the bias
  amplification score with a bootstrap interval; `tmpaud()` the same by
  city and period with each city's range. `calaud()` rolls records up to
  areas and asks whether the order of predicted risk matches the order of
  recorded outcomes (Spearman's rho with the AS 89 and large-sample
  p-values) and whether the mean rank gap of each group's areas departs
  from chance (a permutation test); `scrdsp()` a score's spread across
  groups with the one-way analysis of variance and Welch's test.
  `lipaud()` is the Lipschitz audit of individual fairness. `nsbnd()`
  bounds a true base rate behind a noisily recorded one for noise rates in
  boxes, and says which groups' order survives every admissible noise;
  `impfpr()` is Chouldechova's identity; `rnkres()` the pairs a noisy
  ranking orders. `oaxbd()` is the Oaxaca-Blinder decomposition of a gap
  at either group's coefficients and the threefold one, with delta-method
  standard errors. `fairsim()` simulates records with a known disparity and
  `nordet()` the noisy-OR chance of detection by patrols. `proof()` names
  the Lean theorems behind `nsbnd()`, `impfpr()`, `rnkres()` and
  `oaxbd()`.
* `relchk()` recomputes a reported odds ratio, risk ratio or risk
  difference from its table under every labelling (each reference group,
  each outcome as the event, groups in rows or in columns) and says which
  readings reproduce the value as rounded: the labelling the report
  states, another one, or none. Zero cells give `Inf` or `NaN`, never a
  correction. `proof("relchk")` names the Lean theorems on how the
  contrasts change under relabelling.
* `catchk()` audits categorical inputs before an analysis: storage type,
  levels and counts, missing values, unused and rare levels, labels that
  differ only by case or white space, 0/1 and logical codings, numbers
  that look like category codes, value labels that disagree with the
  values, and the reference level a default would take and why. It can
  check a stated reference level, and changes nothing.
* Disparity tests: `vodtst()`, `orrisk()`, `outtst()`, `dspbm()`,
  `dyadrt()`.
* Forensic statistics: `gtmp()`, `lrpost()`, `cllr()`, `cllrmin()`,
  `lrcal()`, `mixlr()`, `kinidx()`, `patidx()`, `dbsrch()`, `tipplt()`,
  `eceplt()`.
* `clrsmp()` splits a change in a clearance rate into rate and mix effects
  and flags Simpson's paradox.
* `expbnd()` bounds rates and disparities when exposure is known only
  through a proxy.
* `misrd()` gives the risk difference when group membership is misrecorded.
* `rtmhs()` measures regression to the mean at selected hot spots.
* `incap()` gives the crime rate and share prevented under incapacitation.
* `incapr()` gives crimes prevented under desistance and replacement.
* `sellab()` bounds the failure rate of a release rule from released
  defendants alone.
* `losbnd()` bounds the mean time in custody from completed and open stays.
* `sntbnd()` bounds the effect of a sentence on a binary outcome without
  assumptions.
* `hzdsel()` shows how period-two hazards differ after selective depletion.
* Networks: `coofg()` builds the co-offending graph of an incident table;
  `ntdeg()`, `ntbtw()`, `ntcls()` and `ntevc()` give degree, betweenness,
  closeness and eigenvector centrality; `ntcore()` cores, `ntcmp()`
  components, `ntlvn()` and `ntmod()` Louvain communities and modularity,
  `ntbrk()` Gould-Fernandez brokerage and `ntchg()` the Cheeger bounds.
  Results are by node label in order of first appearance.
* Spatial: `hvdist()` and `eqproj()`/`eqinv()` for longitude and latitude;
  `cgstat()`, `cgmedn()`, `webpt()`, `sdell()` and `cvhull()` for centres,
  spread and hulls; `nnidx()`, `ripk()` and `ripenv()` for nearest
  neighbours and K functions; `fzmode()`, `kmpart()`, `nnhclu()` and
  `stacel()` for hot spots; `mantst()` and `stmavg()` in space and time;
  `jtdecay()`, `jtsurf()`, `jtcal()`, `jtacc()` and `bjtc()` for journey to
  crime. Points come as vectors, tables or sf points; longitude and
  latitude only with `lonlat = TRUE`.
* Spatial weights and autocorrelation: `spwt()` (k nearest, distance band,
  contiguity; binary or inverse distance), `gearc()` (Geary's C), `gorg()`
  (the Getis-Ord G), `lmorn()` (local Moran with conditional permutations),
  `locgi()` (local G and G*) and `dirmn()` (directional mean).
* `hdbang()` (head-bang smoothing, weighted or not), `crwalk()` (a
  correlated random walk of a series: fit, expected displacement and
  simulation) and `rnnhcl()` (nearest-neighbour clustering against a
  population at risk). `nnidx()`, `ripk()`, `fzmode()` and `mantst()` take
  `metric = "manhattan"`; `fzmode()`, `kmpart()`, `nnhclu()` and `stacel()`
  take zones as counts at centroids (`count`).
* `fair()` (fairness of a risk tool by group), `tlvc()` and `tlvlr()`
  (trace evidence), `aorwt()` (aoristic weights), `bdig()` and `bnft()`
  (Benford screens), `skp()` (Scott-Knott grouping), `rptc()` (constant
  reporters), `rswp()` and `rshl()` (Rosenbaum bounds), `zgll()`, `zgem()`,
  `mvtll()` and `mvtem()` (group-based trajectories), `tfrp()` (glass
  transfer), `lwss()` (lowess), `chmb()` (tool marks), `gstr()` (grid Gi*)
  and `sptst()` (Andresen's point pattern test).
* `proof()` names the Lean theorem behind each function; the proofs ship
  in `inst/lean`.
* A vector longer than the largest integer is refused with an error, by the
  exported functions and by the compiled entry points alike; the compiled
  code used to truncate its length.
* `eceplt()` gives the calibrated curve its limits at extreme prior odds:
  a calibrated ratio of 0 or infinity used to give `NaN`, or an error.
* `incapr()` refuses a `t0 + s` beyond the length of `lam` however large
  the two are; past the largest integer it used to warn of an integer
  overflow and stop with an unrelated error.
* `sntbnd()` orders the two values of `z` by Unicode code point, so the
  default `treat` and the order of the terms are the same on every machine
  and in every locale (`"B"` sorts before `"a"`).
* Every error message starts with the name of the function that raised it
  (``"outtst: `q` must be at most 1."``); the words after it are unchanged.
* A matrix with more than one row and column is refused where a vector is
  expected, naming the argument; it used to be read as one long vector.
* `eceplt()` labels its vertical axis in bits.
* `proof()` names the theorems now proved for `gtmp()`, `cllr()`,
  `dbsrch()`, `mixlr()`, `kinidx()`, `patidx()`, `misrd()`, `incap()`,
  `lrcal()` and `cllrmin()`. The help pages of `misrd()`, `cllr()` and
  `cllrmin()` say what is proved and what is not.
* `vignette("crimlab")` walks through one example per family.
* Arguments that hold one value per group are matched by name. Named
  inputs whose names come in another order were read by position, a wrong
  answer: `vodtst(c(white = 600, black = 400), c(black = 300, white = 700))`
  read 300 as the white night count. Now every later argument is put in the
  order of the first by name; names that differ, or a name given twice, are
  an error that lists them. Unnamed inputs are still matched by position.
  This holds in `vodtst()`, `dspbm()`, `dyadrt()` (the row and column names
  of `count` too), `clrsmp()`, `expbnd()`, `incap()`, `gtmp()`, `hzdsel()`,
  and for loci in `kinidx()` and `patidx()`; per-unit arguments of
  `outtst()`, `rtmhs()`, `lrpost()`, `sellab()` and `sntbnd()` are matched
  by name when every one has unique names.
* Per-group arguments also accept a list of single numbers and a one-row
  data frame; a one-row or one-column matrix gives its names.
* `vodtst()` takes one table of groups by day and night (a matrix, 2-d
  table or data frame with columns named `day` and `night`) with `night`
  left out.
* `vodtst()`, `dspbm()` and `expbnd()` return a `ref` column marking the
  reference group; a factor `ref` is read by its label.
* `sntbnd()` no longer picks the treatment by sorting the values of `z`:
  `treat` is required unless `z` is logical (`TRUE` is the treatment) or 0
  and 1 (1 is).
* A labelled vector (with a `labels` attribute, as haven makes) is refused
  with a request to decode it first: its numbers are category codes.
* Learning: `enpath()`, the elastic net path by coordinate descent with
  cross-validation over given folds (the objective and scale of glmnet);
  `lgtreg()`, logistic regression; `rfor()`, random forests for
  regression and probabilities.
* Causal estimates: `dmlx()`, double machine learning by cross-fitting
  for the partially linear, interactive, IV, difference-in-differences
  and sample-selection models, with AIPW and TMLE from the same nuisances.
* `enpath()` screens binomial paths by the sequential strong rule (`scr`);
  `rfor()` grows honest trees (`hon`), limits split imbalance (`imb`) and
  takes `mtry` by the rule of grf by default; `dmlx()` adds two-way
  clusters, the experimental DiD scores (`score = "exp"`) and passes `tol`,
  `scr` and `hon` to its learners.
* Policy evaluation on unit-by-period panels. `didfe()` is the two-way
  fixed-effects difference in differences by least squares (alternating
  projections, weights, one-way or two-way cluster-robust errors with the
  small-sample factor of 'fixest'), or by Poisson or negative binomial
  regression for counts, with the estimate leaving out each treated unit,
  the Goodman-Bacon decomposition, and matching on the outcomes before
  treatment. `evtst()` is the event study with the joint test of the leads
  and the raw trends of each adoption cohort. `stgdid()` gives the
  group-time effects of Callaway and Sant'Anna, with the doubly robust,
  IPW and regression estimators for covariates, Wooldridge's regression on
  cohort-period cells, and de Chaisemartin and D'Haultfoeuille's estimator
  for treatments that switch on and off, with cluster-bootstrap errors.
  `synthc()` is the synthetic control of one treated unit with the placebo
  test, and `sdidx()` the synthetic difference in differences with placebo
  or jackknife errors; both solve their weights exactly (Wolfe's
  minimum-norm point). `evtst()` bounds the effects under departures from
  parallel trends (Rambachan and Roth's relative magnitudes) with a
  breakdown value. `wcbx()` gives the restricted wild cluster bootstrap
  (every sign pattern when few), permutation of adoption dates across
  clusters, and simulated power with the minimum detectable effect;
  `wystep()` the Westfall-Young stepdown p-values.
* Spatial crime analysis. `knxnr()` is the Knox test of near repeats: pairs
  by band of distance and time against their exact expectation under
  exchangeable times, with permutation and Poisson p-values. `kdsurf()`
  gives kernel intensity surfaces (the bandwidth and values of
  `MASS::kde2d()`, the edge corrections of `spatstat.explore::density.ppp()`)
  with the smallest areas holding each share of the incidents; `dbscn()`
  DBSCAN clusters as `dbscan::dbscan()` numbers them; `hsper()` hot grid
  cells per period, their overlap and how long each stays hot; `ebshr()`
  the empirical-Bayes shrinkage of hot-spot counts and the fall to expect
  without any programme; `kscan()` Kulldorff's space-time permutation scan
  with secondary clusters. `gmoran()` and `bvmor()` give global and
  bivariate Moran's I with their inference (`spdep::moran.test()`,
  `spdep::moran_bv()`), one variable or a sweep of columns; `sarml()` fits
  the spatial lag, error and combined models by maximum likelihood, and the
  lag model by two-stage least squares, with Durbin terms and the OLS fit
  beside (`spatialreg::lagsarlm()`, `errorsarlm()`, `sacsarlm()`,
  `stsls()`); `spill()` gives direct, indirect and total effects with
  delta-method standard errors (`spatialreg::impacts()`). `crossk()` gives
  the cross-type K function and the kernel segregation test, `qdrtst()` the
  quadrat test, `morsta()` Morisita's index, `dissix()` the dissimilarity,
  entropy, Atkinson and Gini indices with Morrill's adjustment. `pinpl()`
  puts points in polygons and reports the population their regions cover,
  `arwt()` moves counts between two sets of areal units, and `bufas()`
  assigns treated, comparison and excluded places around sites.
  `proof("ebshr")`, `proof("kscan")`, `proof("dissix")`, `proof("spill")`
  and `proof("bvmor")` name the Lean theorems in `inst/lean/Hotspot.lean`.
