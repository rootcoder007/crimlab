# Changelog

## crimlab 0.0.0.9000

- Fixed:
  [`xcser()`](https://rootcoder007.github.io/crimlab/reference/xcser.md)
  used to report r = 0 with an interval, p = 1 and nonzero
  cross-correlations for a constant series whose mean rounds
  (`rep(1, 10)`, `rep(0.1, 10)`); a correlation with a constant series
  is undefined, and
  [`xcser()`](https://rootcoder007.github.io/crimlab/reference/xcser.md)
  and the command line now refuse it, naming the series.

- Fixed: the other statistics that divide by the spread of the data used
  to return numbers for a constant input whose mean rounds
  (`rep(0.1, n)`, `rep(1 / 3, n)`), where they refused a constant with
  an exact mean:
  [`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md),
  [`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md),
  [`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md)
  and
  [`bvmor()`](https://rootcoder007.github.io/crimlab/reference/bvmor.md)
  (Geary’s C, local Moran, G and the bivariate Moran are undefined),
  [`txcal()`](https://rootcoder007.github.io/crimlab/reference/txcal.md)
  (no calibration line),
  [`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)
  (no Imbens-Kalyanaraman bandwidth),
  [`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md)
  (a pair with a constant unit, now skipped) and the Jarque-Bera test of
  [`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md).
  They now refuse it, or skip it, as for any constant.
  `locgi(star = FALSE)` gives `NaN` for the one unit whose value differs
  from all the others (its reference values are constant).
  [`cmpsc()`](https://rootcoder007.github.io/crimlab/reference/cmpsc.md)
  gives no correlation to a stretch that does not vary, as
  [`cor()`](https://rdrr.io/r/stats/cor.html) does, instead of one made
  of the rounding of its running sums.
  [`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md)
  and
  [`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md)
  refuse a covariate column that does not vary (it cannot be told from
  the intercept), naming it, and
  [`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md)
  an outcome that does not vary.

- Fixed:
  [`hsacc()`](https://rootcoder007.github.io/crimlab/reference/hsacc.md)
  used to break ties in the predicted score at the cut-off of `top` by
  the order of the rows, so that reordering the cells changed the hit
  rate, PAI, PEI and RRI (0.571 or 0.536 at `top = 0.2` in one example);
  and `nmax` broke ties in density the same way. The new argument `ties`
  says which tied cells are hot: `"in"` (all, the default), `"out"`
  (none) or `"frac"` (each by the share of their area still needed);
  `nmax` takes the larger of equally dense cells first. The result no
  longer depends on the order of the rows; without ties at the cut-off
  it is as before.

- Fixed: the point patterns and distance functions used to run without
  end, or until the machine was stopped, on extreme radii and
  magnitudes: a radius or distance of 5e-324
  ([`ripk()`](https://rootcoder007.github.io/crimlab/reference/ripk.md),
  [`ripenv()`](https://rootcoder007.github.io/crimlab/reference/ripenv.md),
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md),
  [`knxnr()`](https://rootcoder007.github.io/crimlab/reference/knxnr.md)
  and the other functions on the shared grid of points), a radius of 1e9
  ([`ripk()`](https://rootcoder007.github.io/crimlab/reference/ripk.md),
  [`ripenv()`](https://rootcoder007.github.io/crimlab/reference/ripenv.md)),
  `opteffort(W = 5e-324)`,
  [`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md)
  and
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  on counts of 1e50 and
  [`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md)
  on coordinates of 1e200. The grid now starts from the smallest normal
  double and its searches stop at its edge;
  [`opteffort()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)
  refuses a `W` whose densities underflow or overflow; the goodness of
  fit of
  [`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md)’s
  discrete fits looks only at the jumps of the data;
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  refuses counts above 2^53 for its exact test;
  [`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md)
  and
  [`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md)
  refuse points whose squared distances overflow; and the optimisers
  stop on an infinite start or search direction. Results for other
  inputs are unchanged.

- Fixed:
  [`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md)
  and
  [`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)
  used to refuse a series that does not vary once differenced (a
  constant, a straight line) in the optimiser’s words (“initial value in
  ‘vmmin’ is not finite”); they now say “`y` does not vary once
  differenced: the model has nothing to fit.”, as does the command line.

- Fixed: the command line’s
  [`fylbl()`](https://rootcoder007.github.io/crimlab/reference/yoylbl.md)
  converted an `endyear` of 1e300, `Inf` or `NaN` to a whole number with
  undefined behaviour; like R it now gives `NA` for a year outside the
  range of an integer.
  [`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md)
  counted `days / dt` steps as an integer before checking them (1e16 at
  `dt = 1e-9`), and the periodic baseline of the Hawkes functions its
  integration steps; both are now counted in double first.

- Fixed: the command line and `clab_jrun()` used to copy text that is
  not valid UTF-8 from an argument or a CSV cell into their output,
  which was then not valid JSON; they now refuse it, naming the
  argument, in the words of R (“`group` holds text that is not valid
  UTF-8.”). A NUL byte in a CSV cell used to end the cell’s number
  silently (`4<NUL>0` was read as 4); it is now refused with its row and
  column.

- Fixed: a function that calls another
  ([`runfx()`](https://rootcoder007.github.io/crimlab/reference/exprun.md)
  calls
  [`exprun()`](https://rootcoder007.github.io/crimlab/reference/exprun.md)
  and an estimator,
  [`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md)
  [`gplume()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
  [`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
  [`crfno2()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md),
  [`rsvsmp()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md)
  [`rsvidx()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md),
  [`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md)
  [`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md),
  [`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md)
  the writers, and the exposure functions
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md))
  used to pass on the inner function’s refusals and warnings under the
  inner name (“gplume: …”); they now carry the name of the function
  called
  ([`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md),
  [`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md),
  [`hetfx()`](https://rootcoder007.github.io/crimlab/reference/hetfx.md),
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md),
  [`lftab()`](https://rootcoder007.github.io/crimlab/reference/lftab.md),
  [`plcday()`](https://rootcoder007.github.io/crimlab/reference/plcday.md),
  [`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md),
  [`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md),
  [`prechk()`](https://rootcoder007.github.io/crimlab/reference/prechk.md),
  [`rsvsmp()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md),
  [`runfx()`](https://rootcoder007.github.io/crimlab/reference/exprun.md),
  [`typefx()`](https://rootcoder007.github.io/crimlab/reference/typefx.md),
  [`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md),
  [`yoyhtml()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md),
  [`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md)).

- Fixed: the choice arguments of the exposure functions
  ([`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
  [`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md),
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md),
  [`expgmm()`](https://rootcoder007.github.io/crimlab/reference/expgmm.md),
  [`hetfx()`](https://rootcoder007.github.io/crimlab/reference/hetfx.md),
  [`invsn()`](https://rootcoder007.github.io/crimlab/reference/invsn.md),
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md),
  [`plcday()`](https://rootcoder007.github.io/crimlab/reference/plcday.md),
  [`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md),
  [`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md),
  [`typefx()`](https://rootcoder007.github.io/crimlab/reference/typefx.md),
  [`updwn()`](https://rootcoder007.github.io/crimlab/reference/updwn.md),
  [`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md))
  used to be refused by
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) in base R’s
  words (“‘arg’ should be one of …”); they are now refused like every
  other choice (“ppmlfe: `se` must be one of …”).

- Fixed:
  [`hwpwr()`](https://rootcoder007.github.io/crimlab/reference/hwpwr.md)
  used to report as measured a RAPL counter that fell by more than a
  wrap could explain (a reset read as 257 kJ) or read above its range,
  and read counter text that is not a whole decimal number (“0x10”,
  “1.5e3”, an embedded NUL, 26 digits); such rows are now not measured,
  with the reason in `why`. A FIFO in place of a counter used to block
  the call for ever; files are now read only when regular, opened
  without blocking. A library that is not NVML, or an `nvidia-smi` that
  fails, no longer adds a GPU row. `dt` is at most 604800 (a week) and
  `pid` at most 2147483647 (`pid = 2^31` used to be accepted with a
  warning).

- Environmental crime, with offences against the environment as the
  outcome: the LMDI decomposition of a count over the factors of a Kaya
  identity, over time and against a reference region
  ([`lmdix()`](https://rootcoder007.github.io/crimlab/reference/lmdix.md),
  after Duan et al. 2026), the two-way fixed-effects regression of a
  count on drivers with an interaction, standardised coefficients,
  variance inflation factors, simple slopes and the Jarque-Bera test of
  the residuals
  ([`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md)),
  and Pesaran’s CD test of cross-sectional dependence
  ([`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md),
  as Andreatta et al. 2023 use it). They are in `clab_jrun()`, the
  command line and Python.

- Dynamic and heterogeneous panels, for the models Andreatta et
  al. (2023) fit to waste crime: difference and system GMM
  ([`gmmpnl()`](https://rootcoder007.github.io/crimlab/reference/gmmpnl.md),
  Arellano and Bond 1991, Arellano and Bover 1995, Blundell and Bond
  1998, with the instruments of Roodman 2009), one-step and two-step,
  first differences or forward orthogonal deviations, GMM-style
  instruments with lag limits and collapsing, IV-style instruments, the
  robust or conventional one-step variance and Windmeijer’s corrected or
  the uncorrected two-step variance, small-sample t and F statistics, a
  joint Wald or F test, the Arellano-Bond tests, the Sargan, Hansen and
  difference-in-Hansen tests and the instrument count; Wooldridge’s test
  of serial correlation
  ([`wsctst()`](https://rootcoder007.github.io/crimlab/reference/wsctst.md),
  Drukker 2003); and the common correlated effects mean group and pooled
  estimators
  ([`ccepnl()`](https://rootcoder007.github.io/crimlab/reference/ccepnl.md),
  Pesaran 2006), on balanced and unbalanced panels, with a choice of the
  averages taken out. They are in `clab_jrun()`, the command line and
  Python.

- Panels with a repeated unit and period are found by sorting the rows,
  not with a grid of units by periods, which for many units and periods
  could ask for more memory than the rows hold.
  [`pcdtst()`](https://rootcoder007.github.io/crimlab/reference/pcdtst.md)
  and
  [`lmdix()`](https://rootcoder007.github.io/crimlab/reference/lmdix.md)
  refuse a panel whose units by periods pass 2^26 cells.

- Spatial epidemiology and surveillance: Bayes outbreak detection
  ([`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md)),
  exposure within a buffer and by kernel
  ([`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
  [`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md)),
  the Leroux and BYM2 area structures
  ([`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
  [`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md)),
  ecological Poisson, negative binomial and zero-inflated regression
  with expected counts
  ([`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md)),
  eigenvector spatial filtering
  ([`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md)),
  the metapopulation SIR model
  ([`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)),
  Weibull frailty and promotion-time cure models with CAR area effects
  ([`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
  [`sfcure()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md))
  and the spatial encoding of graph attention
  ([`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md)).
  They are in `clab_jrun()`, the command line and Python.

- Geostatistics on planar coordinates: the empirical semivariogram
  ([`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md)),
  variogram models fitted by maximum likelihood
  ([`vgfit()`](https://rootcoder007.github.io/crimlab/reference/vgfit.md)),
  ordinary kriging
  ([`okrge()`](https://rootcoder007.github.io/crimlab/reference/okrge.md)),
  block kriging
  ([`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md))
  and its leave-one-out cross-validation
  ([`krigcv()`](https://rootcoder007.github.io/crimlab/reference/krigcv.md)),
  polynomial trend surfaces
  ([`trsurf()`](https://rootcoder007.github.io/crimlab/reference/trsurf.md)),
  the spatial linear model by profile likelihood
  ([`sglmx()`](https://rootcoder007.github.io/crimlab/reference/sglmx.md))
  and the spatial linear mixed model by REML
  ([`smixm()`](https://rootcoder007.github.io/crimlab/reference/smixm.md)),
  with a spatial data object
  ([`spobj()`](https://rootcoder007.github.io/crimlab/reference/spobj.md),
  [`spxy()`](https://rootcoder007.github.io/crimlab/reference/spxy.md),
  [`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md)),
  the specification of a web map
  ([`lfspec()`](https://rootcoder007.github.io/crimlab/reference/lfspec.md))
  and density-weighted sampling of points
  ([`spsamp()`](https://rootcoder007.github.io/crimlab/reference/spsamp.md)).
  They are in `clab_jrun()`, the command line and Python.

- The energy and carbon footprint of a computation
  ([`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md)),
  modelled from running time, utilisation and memory (Green Algorithms)
  or measured from the CPU’s RAPL counters
  ([`raplok()`](https://rootcoder007.github.io/crimlab/reference/raplok.md)),
  with the measured energy of other devices added from a table; the
  carbon intensity of electricity by location, offline
  ([`co2int()`](https://rootcoder007.github.io/crimlab/reference/co2int.md),
  [`co2inttab()`](https://rootcoder007.github.io/crimlab/reference/co2int.md)),
  the location found without network access
  ([`dtctloc()`](https://rootcoder007.github.io/crimlab/reference/dtctloc.md))
  and everyday equivalents
  ([`fpeqv()`](https://rootcoder007.github.io/crimlab/reference/fpeqv.md)).
  All but
  [`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md),
  which evaluates an R expression, are in `clab_jrun()`, the command
  line and Python.

- Air pollution, for environmental criminology: the exposure-response
  estimate with a bootstrap sensitivity interval
  ([`exprspplr()`](https://rootcoder007.github.io/crimlab/reference/exprspplr.md)),
  by the cross-fitted partially linear regression in the compiled core
  or any estimator with its interface, and the whole pollution-to-health
  pipeline from demo data, a (compressed) CSV or scalars to a report and
  its text
  ([`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md),
  [`polrpttxt()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)).
  [`exprspplr()`](https://rootcoder007.github.io/crimlab/reference/exprspplr.md)
  with its built-in estimator is in `clab_jrun()`, the command line and
  Python; the other two read files and write text in R.

- Scaling and ideal points, for judges, courts, legislators and survey
  respondents: vote agreement
  ([`votagr()`](https://rootcoder007.github.io/crimlab/reference/votagr.md)),
  Aldrich-McKelvey scaling
  ([`amscal()`](https://rootcoder007.github.io/crimlab/reference/amscal.md))
  and its Bayesian form
  ([`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md)),
  Basic Space scaling
  ([`bbscal()`](https://rootcoder007.github.io/crimlab/reference/bbscal.md)),
  classical, SMACOF and nonmetric MDS
  ([`cmdsf()`](https://rootcoder007.github.io/crimlab/reference/cmdsf.md),
  [`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md),
  [`nmdsf()`](https://rootcoder007.github.io/crimlab/reference/nmdsf.md),
  [`dblctr()`](https://rootcoder007.github.io/crimlab/reference/dblctr.md),
  [`mdsfit()`](https://rootcoder007.github.io/crimlab/reference/mdsfit.md)),
  Bayesian MDS and unfolding
  ([`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md),
  [`bunfld()`](https://rootcoder007.github.io/crimlab/reference/bunfld.md)),
  unfolding by least squares and by SMACOF
  ([`mlsmu()`](https://rootcoder007.github.io/crimlab/reference/mlsmu.md),
  [`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md),
  [`unfstr()`](https://rootcoder007.github.io/crimlab/reference/unfstr.md),
  [`ipsrec()`](https://rootcoder007.github.io/crimlab/reference/ipsrec.md)),
  INDSCAL
  ([`indscl()`](https://rootcoder007.github.io/crimlab/reference/indscl.md)),
  Procrustes rotation
  ([`procrs()`](https://rootcoder007.github.io/crimlab/reference/procrs.md)),
  NOMINATE utilities, likelihood, fit and bootstrap
  ([`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md),
  [`nompr()`](https://rootcoder007.github.io/crimlab/reference/nompr.md),
  [`nomll()`](https://rootcoder007.github.io/crimlab/reference/nomll.md),
  [`dwnom()`](https://rootcoder007.github.io/crimlab/reference/dwnom.md),
  [`nombt()`](https://rootcoder007.github.io/crimlab/reference/nombt.md)),
  alpha-NOMINATE
  ([`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md)),
  Optimal Classification and its ordered form
  ([`optcls()`](https://rootcoder007.github.io/crimlab/reference/optcls.md),
  [`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md)),
  binary IRT by EM
  ([`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md))
  and by Gibbs sampling
  ([`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md),
  [`irtlik()`](https://rootcoder007.github.io/crimlab/reference/irtlik.md),
  [`irtpst()`](https://rootcoder007.github.io/crimlab/reference/irtpst.md)),
  ordinal and dynamic IRT
  ([`ordirt()`](https://rootcoder007.github.io/crimlab/reference/ordirt.md),
  [`dynirt()`](https://rootcoder007.github.io/crimlab/reference/dynirt.md)),
  anchoring vignettes
  ([`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md)),
  the bootstrap of scaling positions
  ([`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md)),
  normal vectors and cutting lines
  ([`nrmvec()`](https://rootcoder007.github.io/crimlab/reference/nrmvec.md),
  [`cutlns()`](https://rootcoder007.github.io/crimlab/reference/cutlns.md))
  and Wordfish
  ([`wrdfsh()`](https://rootcoder007.github.io/crimlab/reference/wrdfsh.md)).
  The samplers draw from R’s generator seeded as
  [`set.seed()`](https://rdrr.io/r/base/Random.html)
  ([`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md)
  from Philox), so a chain is the same in every front end. All are in
  `clab_jrun()`, the command line and Python.

- Counts over time and across categories: ARIMA and seasonal ARIMA
  forecasts by exact maximum likelihood
  ([`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
  [`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)),
  the Ornstein-Uhlenbeck fit, simulation and Fokker-Planck density of a
  count series
  ([`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
  [`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md)),
  Lotka-Volterra rates of a yearly series
  ([`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md))
  and the correlation of counts per unit between categories
  ([`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md)).
  They are in `clab_jrun()`, the command line and Python.

- Regression discontinuity: local-linear sharp and fuzzy estimates with
  the bandwidth of Imbens and Kalyanaraman (2012) or one given, three
  kernels, conventional and robust intervals (Calonico, Cattaneo and
  Titiunik 2014, with the pilot bandwidth equal to the main one),
  heteroskedasticity-robust or clustered errors and placebo cutoffs
  ([`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md));
  the density test of McCrary (2008)
  ([`rddens()`](https://rootcoder007.github.io/crimlab/reference/rddens.md)).

- Exposure and crime, checked against the papers: the Newey-West lag of
  `se = "nw"` when `lag` is not given is now the automatic bandwidth of
  Newey and West (1994, equation 2.2) instead of their lag-selection
  parameter `n`, so results with that default change; two-way clustering
  gains the factor of each component (`twc = "each"`) and the eigenvalue
  fix (`psd = "eig"`) of Cameron, Gelbach and Miller (2011) in
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md)
  and
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md).
  The first-stage table of
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  gains the Kleibergen-Paap rk statistic for several endogenous columns
  (it was `NA`) and the effective degrees of freedom and 5% critical
  value of the simplified test of Montiel Olea and Pflueger (2013); the
  effective F itself was checked against their equation 3 and is
  unchanged.

- [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md)
  and
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  absorb any number of sets of fixed effects (a third set such as hour
  or region-by-year, which several of the papers use); one or two sets
  are computed exactly as before.

- [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  gains single-day lags (`lagv`, `lag`) and natural-spline controls
  (`ns`), the specification of the nationwide case-crossover study of
  Han et al. (2026), with a worked example of its main table.

- `smkiv(cf = "childs")`: smoke PM2.5 as Childs et al. (2022) define it,
  PM2.5 above the three-year median of the same place and month on
  smoke-free days, on smoke days.

- `ppmlfe(se = "quasi")`, the quasi-Poisson variance with the Pearson
  dispersion, the overdispersion adjustment Armstrong, Gasparrini and
  Tobias (2014) give for the case-crossover design of
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md).

- [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  reports the Sargan statistic and the cluster-robust Hansen J of the
  two-step GMM estimate when there are more instruments than endogenous
  columns (`overid`, `gmm`), the cross-definition over-identification
  test the review of Yamada (2026) recommends.

- `ppmlfe(sep = "sst")`: the check of Santos Silva and Tenreyro (2010)
  for Poisson estimates that do not exist, dropping the regressors that
  are collinear on the rows with positive counts. The iterative
  rectifier of Correia, Guimaraes and Zylkin, which drops rows instead,
  is not built: their 2020 paper names it without stating it.

- Exposure and crime, a count with an endogenous exposure by the
  exponential-mean generalised method of moments of Mullahy (1997) and
  Windmeijer and Santos Silva (1997): multiplicative or additive
  residual, one-step, two-step or iterated weighting, cluster-robust
  errors and Hansen’s J
  ([`expgmm()`](https://rootcoder007.github.io/crimlab/reference/expgmm.md)).

- Exposure and crime, the designs of 2019-2026: smoke instruments from
  plume data, raw or adjusted to the plumes that reached the ground
  ([`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md));
  runs of exposed days and the effect by day in the run
  ([`exprun()`](https://rootcoder007.github.io/crimlab/reference/exprun.md),
  [`runfx()`](https://rootcoder007.github.io/crimlab/reference/exprun.md));
  effects by a moderator with the test of equality
  ([`hetfx()`](https://rootcoder007.github.io/crimlab/reference/hetfx.md));
  the time-stratified case-crossover design
  ([`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md));
  two exposures together as a surface
  ([`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md));
  exposure in bins with tests of non-linearity
  ([`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md));
  and effects by offence type with a multiple-testing correction
  ([`typefx()`](https://rootcoder007.github.io/crimlab/reference/typefx.md)).
  Every choice the papers leave open is an option of the estimators,
  with the standard errors of each kind, the three first-stage
  statistics and the plug-in of the count instrument.

- Exposure and crime, for counts by place and day against a daily
  exposure: Poisson pseudo-maximum likelihood with one or two sets of
  absorbed fixed effects, an offset and one- or two-way clustered errors
  ([`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md));
  two-stage least squares with the same effects and the first-stage Wald
  statistic
  ([`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md));
  the Poisson control function with a cluster bootstrap
  ([`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md));
  wind-direction and inversion instruments
  ([`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md),
  [`invsn()`](https://rootcoder007.github.io/crimlab/reference/invsn.md));
  the irrelevant-day placebo
  ([`plcday()`](https://rootcoder007.github.io/crimlab/reference/plcday.md));
  the upwind and downwind design around a road or a point source
  ([`updwn()`](https://rootcoder007.github.io/crimlab/reference/updwn.md));
  and the join of exposure to counts
  ([`expjn()`](https://rootcoder007.github.io/crimlab/reference/expjn.md)).
  The methods are those of Bondy, Roth and Sager (2020) and of
  Herrnstadt, Heyes, Muehlegger and Saberian (2019, 2021).

- Air pollution, for environmental criminology: Pasquill-Gifford
  dispersion coefficients
  ([`pgsig()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md)),
  Briggs plume rise
  ([`briggs()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md)),
  the Gaussian plume and puff
  ([`gplume()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
  [`gspuff()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
  the field
  [`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md)),
  an advection-diffusion grid
  ([`advdif2d()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md))
  and a particle walk
  ([`lagrprt()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md));
  the concentration-response functions of PM2.5 and NO2
  ([`crfpm25()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md),
  [`crfno2()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md)),
  the attributable fraction
  ([`attrfrac()`](https://rootcoder007.github.io/crimlab/reference/attrfrac.md)),
  deaths avoided
  ([`mortdisp()`](https://rootcoder007.github.io/crimlab/reference/mortdisp.md)),
  the burden of a population and of areas
  ([`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md),
  [`polbrdarea()`](https://rootcoder007.github.io/crimlab/reference/polbrdarea.md))
  and the concentration index of exposure by income
  ([`expconcidx()`](https://rootcoder007.github.io/crimlab/reference/expconcidx.md)).

- Describing a table before trusting it: inline histograms
  ([`inlhist()`](https://rootcoder007.github.io/crimlab/reference/inlhist.md)),
  frequency tables
  ([`freqtbl()`](https://rootcoder007.github.io/crimlab/reference/freqtbl.md)),
  pairwise correlations
  ([`cortbl()`](https://rootcoder007.github.io/crimlab/reference/cortbl.md),
  [`topcors()`](https://rootcoder007.github.io/crimlab/reference/topcors.md)),
  empty and constant rows and columns
  ([`drpempty()`](https://rootcoder007.github.io/crimlab/reference/drp.md),
  [`drpconst()`](https://rootcoder007.github.io/crimlab/reference/drp.md)),
  clean column names
  ([`clncolnms()`](https://rootcoder007.github.io/crimlab/reference/clncolnms.md)),
  multivariate outliers
  ([`mahaoutl()`](https://rootcoder007.github.io/crimlab/reference/mahaoutl.md)),
  the shape of missingness
  ([`missruns()`](https://rootcoder007.github.io/crimlab/reference/missruns.md),
  [`missgrid()`](https://rootcoder007.github.io/crimlab/reference/missgrid.md),
  [`misspat()`](https://rootcoder007.github.io/crimlab/reference/misspat.md),
  [`misssmry()`](https://rootcoder007.github.io/crimlab/reference/misssmry.md))
  and Little’s test of missing completely at random
  ([`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md)),
  duplicated rows
  ([`duprows()`](https://rootcoder007.github.io/crimlab/reference/duprows.md)),
  a column report
  ([`profcols()`](https://rootcoder007.github.io/crimlab/reference/profcols.md)),
  candidate distributions fitted by maximum likelihood and ranked by AIC
  ([`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md))
  and Cochran’s Q for distinct groups before pooling
  ([`dstntst()`](https://rootcoder007.github.io/crimlab/reference/dstntst.md)).
  They are in `clab_jrun()`, the command line and Python.

- One-pass summaries of data larger than memory: moment accumulators
  ([`onlsmry()`](https://rootcoder007.github.io/crimlab/reference/onlsmry.md),
  [`smryupd()`](https://rootcoder007.github.io/crimlab/reference/onlsmry.md),
  [`smrymrg()`](https://rootcoder007.github.io/crimlab/reference/onlsmry.md),
  [`smrystat()`](https://rootcoder007.github.io/crimlab/reference/onlsmry.md)),
  reservoir sampling
  ([`rsvidx()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md),
  [`rsvsmp()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md))
  and the HyperLogLog distinct count
  ([`dstnskch()`](https://rootcoder007.github.io/crimlab/reference/dstnskch.md),
  [`dstncnt()`](https://rootcoder007.github.io/crimlab/reference/dstnskch.md),
  [`skchmrg()`](https://rootcoder007.github.io/crimlab/reference/dstnskch.md)).

- Search theory for search and rescue planning: lateral range curves
  ([`ltrlrng()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md))
  and sweep width
  ([`swpwdth()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)),
  the probability of detection from coverage
  ([`srchpod()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md))
  and of a parallel sweep with navigation error
  ([`pswppod()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)),
  the probability of success with Bayesian updating of containment
  ([`srchscs()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)),
  datum probability maps
  ([`datumpmap()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)),
  the optimal allocation of effort
  ([`opteffort()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md),
  [`optcrcsrch()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)),
  Frost’s optimal search squares
  ([`optsqr()`](https://rootcoder007.github.io/crimlab/reference/optsqr.md))
  and the list of his book’s figures
  ([`srchfig()`](https://rootcoder007.github.io/crimlab/reference/srchfig.md)).
  They are in `clab_jrun()`, the command line and Python.

- First Nations, Metis and Inuit health and wellbeing, a core part of
  the package beside the justice system it concerns: Statistics Canada’s
  published tables read with their quality flags
  ([`sctidy()`](https://rootcoder007.github.io/crimlab/reference/sctidy.md),
  [`scests()`](https://rootcoder007.github.io/crimlab/reference/scests.md),
  [`qltyrcd()`](https://rootcoder007.github.io/crimlab/reference/qltyrcd.md),
  [`limfrom()`](https://rootcoder007.github.io/crimlab/reference/limfrom.md)),
  compared group by group and as whole rated scales
  ([`indgap()`](https://rootcoder007.github.io/crimlab/reference/indgap.md),
  [`ordcmp()`](https://rootcoder007.github.io/crimlab/reference/ordcmp.md)),
  standardised by age
  ([`stdprev()`](https://rootcoder007.github.io/crimlab/reference/stdprev.md)),
  wellness beside distress
  ([`welprof()`](https://rootcoder007.github.io/crimlab/reference/welprof.md),
  [`dualcont()`](https://rootcoder007.github.io/crimlab/reference/dualcont.md),
  [`holbal()`](https://rootcoder007.github.io/crimlab/reference/holbal.md),
  [`frmaud()`](https://rootcoder007.github.io/crimlab/reference/frmaud.md)),
  the First Nations Mental Wellness Continuum
  ([`mhcont()`](https://rootcoder007.github.io/crimlab/reference/mhcont.md)),
  healing as change
  ([`healchg()`](https://rootcoder007.github.io/crimlab/reference/healchg.md)),
  connection across places
  ([`connassoc()`](https://rootcoder007.github.io/crimlab/reference/connassoc.md)),
  community-weighted indices
  ([`commidx()`](https://rootcoder007.github.io/crimlab/reference/commidx.md))
  and the Community Well-Being index
  ([`cwbscr()`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md));
  the generations of residential school attendance and their
  dose-response
  ([`genexpo()`](https://rootcoder007.github.io/crimlab/reference/genexpo.md),
  [`exptrnd()`](https://rootcoder007.github.io/crimlab/reference/exptrnd.md));
  the 2021 Census Indigenous Population Profile
  ([`cnsprof()`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md),
  [`profshr()`](https://rootcoder007.github.io/crimlab/reference/profshr.md),
  [`cnsstd()`](https://rootcoder007.github.io/crimlab/reference/cnsstd.md),
  [`cnscwb()`](https://rootcoder007.github.io/crimlab/reference/cnscwb.md));
  and the Calls to Action
  ([`trcind()`](https://rootcoder007.github.io/crimlab/reference/trcind.md),
  [`trchlthgap()`](https://rootcoder007.github.io/crimlab/reference/trchlthgap.md),
  [`trccmpl()`](https://rootcoder007.github.io/crimlab/reference/trccmpl.md),
  [`trcprog()`](https://rootcoder007.github.io/crimlab/reference/trcprog.md)).
  They are in `clab_jrun()`, the command line and Python.

- [`polex()`](https://rootcoder007.github.io/crimlab/reference/polex.md)
  makes the simulated policy panel the examples of the policy functions
  use: twenty districts, 36 months and a staged reform with a known
  effect. The same seed gives the same panel in every front end.

- [`thrtst()`](https://rootcoder007.github.io/crimlab/reference/thrtst.md):
  when the quasi-Newton steps of the maximum-likelihood fit end short of
  the gradient criterion, up to 20 Newton steps now follow, so the fit
  converges on flat likelihoods on every architecture (arm64 builds
  reported the self-test’s fit as not converged). Fits that met the
  criterion give the same results as before.

- [`thrtst()`](https://rootcoder007.github.io/crimlab/reference/thrtst.md)
  gains `fit = "bayes"`: the threshold test fitted by Bayesian
  inference, as Simoiu, Corbett-Davies and Goel (2017) and Pierson,
  Corbett-Davies and Goel (2017) fit it, with a threshold for every
  group and area and the hierarchical priors that Simoiu et al. state,
  on the discriminant distributions. The posterior is sampled in C++ by
  the No-U-Turn sampler (`chains`, `iter`, `seed`), starting near the
  maximum-likelihood fit. `thr` gives the posterior means and intervals
  of the thresholds and of their differences from the reference group,
  and `post` the summaries, effective sample size and R-hat of every
  parameter. Chains with an R-hat of 1.05 or more, or fewer than 100
  effective draws of a threshold, are refused. The fit agrees with
  Stan’s sampler on the same model. `fit = "ml"` stays the default. The
  two fits answer under different models and should be reported as such,
  not as checks of one another: with a threshold for every group and
  area under a hierarchical prior, on data made with one threshold per
  group the Bayesian fit can report thresholds far from the
  maximum-likelihood ones and from the generating values (three groups
  in ten areas, 20,000 stops a cell: about 0.13, 0.07 and 0.14 against
  0.27, 0.20 and 0.31), which Stan reproduces on the same model.

- [`thrtst()`](https://rootcoder007.github.io/crimlab/reference/thrtst.md)
  gains `rank = "joint"`, the rank-order test of Anwar and Fang (2006)
  as the paper states it: the ranking of all officer groups is tested at
  once by their one-sided z statistics, and the null of no relative
  prejudice is rejected when the significant orders form a cycle; the
  common ranking is given otherwise. The pairwise form stays the
  default. With `off`, the new element `mono` holds their test of
  monolithic behaviour (Pearson’s chi-square of equal rates across
  officer groups, by group of the people stopped and measure).

- [`ibench()`](https://rootcoder007.github.io/crimlab/reference/ibench.md)
  gains the two methods of Ridgeway and MacDonald (2009) it left out:
  `dr = TRUE` gives the doubly robust benchmark (the propensity-weighted
  logistic regression of the outcome on the officer’s indicator and the
  circumstances, with a sandwich standard error), and `local = TRUE`
  Efron’s local false discovery rate of each officer, with the empirical
  null fitted as the package locfdr does by default (column `lfdr`,
  element `null`), flagging by it. Benjamini and Hochberg’s `q` stays.
  The defaults give the results as before. The tests check the doubly
  robust benchmark against [`glm()`](https://rdrr.io/r/stats/glm.html)
  and a sandwich by hand, and the local rate against locfdr (reference
  values from `dev/ref/lfdr.R`).

- [`lnsize()`](https://rootcoder007.github.io/crimlab/reference/lnsize.md)
  follows Malpass (1981) for `esize`: members nobody chose (null foils)
  leave the nominal size and the chance expectation, as Tredoux (1996,
  Table 4.2) shows. Before, `esize` did not adjust; that value is now
  `esizeb`. The two differ only when a member was never chosen. The
  interval of Tredoux’s `E'` (`etrlo`, `etrhi`) now takes the ends of
  the interval of `I = 1 - sum p^2` to `1 / (1 - I)`, as Tredoux does,
  in place of a symmetric delta-method interval around `E'`; the ends
  move for every lineup (for 40 mock witnesses choosing 15, 5, 6, 4, 7
  and 3 of six, E’ 4.44: from 2.98 to 5.90 before, 3.35 to 6.62 now). A
  choice of 0 records a mock witness who chose no one: it counts for
  `psusp` and `fsize` (as in Wells, Leippe and Ostrom 1979) and not for
  the spread of choices; `nm` is the number who chose a member. The
  tests reproduce the worked examples of Tredoux (1996, 1999).

- [`lnroc()`](https://rootcoder007.github.io/crimlab/reference/lnroc.md)
  continues a curve that ends before `far` by the straight line to the
  point where every lineup counts as an identification, as pROC does, in
  place of a flat line; this changes `pauc` only when `far` is given
  above a procedure’s largest false alarm rate (the default never is).
  The tests reproduce Table 3 and the partial areas of Mickes, Flowe and
  Wixted (2012).

- [`wdqx()`](https://rootcoder007.github.io/crimlab/reference/wdqx.md)
  now names the crimes prevented in the treated area and the buffer
  `tne`, the total net effect of Guerette (2009, Appendix B, equation
  4), in place of `net` (same value), and adds the gross effect `ge` and
  the net effect `ne` of the same appendix. Its help page gives the
  formulas and Guerette’s worked example; the tests reproduce it.

- Every result that holds something to draw has a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method on
  base graphics: estimates with intervals as dot-and-interval charts,
  effects by time since adoption, paths, survival and ROC curves, Lorenz
  curves, bars, matrices with their numbers in the cells, and surfaces.
  Each returns, invisibly, the data it drew.
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md)
  gives the colour-blind-safe palettes they share,
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md)
  their session look and
  [`yoypals()`](https://rootcoder007.github.io/crimlab/reference/blpal.md)
  the palettes of change tables; series differ in line type or point
  shape as well as colour. The general plots come with them: Lorenz and
  funnel plots
  ([`plotlrnz()`](https://rootcoder007.github.io/crimlab/reference/plotlrnz.md),
  [`plotfnl()`](https://rootcoder007.github.io/crimlab/reference/plotlrnz.md)),
  distributions by group
  ([`plotrdg()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotcdfs()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotvln()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plothex()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotgain()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)),
  series over time
  ([`plotcal()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plothrz()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plotmult()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plotbump()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md))
  and event studies in a grid
  ([`plotevgrid()`](https://rootcoder007.github.io/crimlab/reference/plotevgrid.md)).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an
  `hwkres` result now returns the quantiles it drew.

- Change tables from
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  (or any table with its columns) can be labelled
  ([`yoylbl()`](https://rootcoder007.github.io/crimlab/reference/yoylbl.md),
  [`fylbl()`](https://rootcoder007.github.io/crimlab/reference/yoylbl.md)),
  summarised
  ([`yoysmry()`](https://rootcoder007.github.io/crimlab/reference/yoysmry.md))
  and written: CSV, TSV and JSON with every digit and every code as it
  is
  ([`yoycsv()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md),
  [`yoytsv()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md),
  [`yoyjs()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md)),
  Markdown
  ([`yoymd()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md)),
  a self-contained HTML page
  ([`yoyhtml()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md))
  and a PDF written by the package itself
  ([`yoypdf()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md)),
  the format picked from the file name by
  [`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md).

- The drift battery: the two-sample Kolmogorov-Smirnov test
  ([`drftks()`](https://rootcoder007.github.io/crimlab/reference/drftks.md)),
  the population stability index and Jensen-Shannon divergence
  ([`drftpsi()`](https://rootcoder007.github.io/crimlab/reference/drftpsi.md))
  and the chi-square tests of a categorical column against a known
  distribution and of two samples’ homogeneity
  ([`drftchisq()`](https://rootcoder007.github.io/crimlab/reference/drftchisq.md),
  [`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md)),
  each result keeping its inputs in the attribute `lakner`. They check a
  new release of a published table against the last one before
  year-over-year change is computed from it. They are in `clab_jrun()`,
  the command line and Python (there without the Monte Carlo branch of
  [`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md)).

- [`rgnmapintg()`](https://rootcoder007.github.io/crimlab/reference/rgnmapintg.md),
  [`rgnmapcmp()`](https://rootcoder007.github.io/crimlab/reference/rgnmapcmp.md)
  and
  [`rgnmapalt()`](https://rootcoder007.github.io/crimlab/reference/rgnmapalt.md)
  check a crosswalk from units to regions: each unit once and with a
  known region, a rebuilt map against the published one, and the map
  against an assignment by another route with the known disagreements
  named. They are in `clab_jrun()`, the command line and Python.

- The categorical-integrity checks: codes decoded by an explicit
  dictionary or by their value labels
  ([`deccodes()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
  [`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md)),
  recodes and relabels by name with every value mapped
  ([`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
  [`relbl()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md)),
  declared levels
  ([`grdlvls()`](https://rootcoder007.github.io/crimlab/reference/grdlvls.md)),
  a 0/1 indicator that refuses a factor
  ([`grdbin()`](https://rootcoder007.github.io/crimlab/reference/grdbin.md)),
  a recode proved by its cross-tabulation
  ([`vrfrcd()`](https://rootcoder007.github.io/crimlab/reference/vrfrcd.md)),
  counts checked against those a release published, with the label
  permutation that would explain a mismatch
  ([`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md)),
  and a variable checked after it crossed from one program to another
  ([`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md)).
  They are in `clab_jrun()`, the command line and Python.

- [`shrs()`](https://rootcoder007.github.io/crimlab/reference/shrs.md)
  gives each count’s share of a total in percent with Wilson’s score
  interval, from a vector or from a table summed by groups. It is in
  `clab_jrun()`, the command line and Python.

- [`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md),
  [`asciifb()`](https://rootcoder007.github.io/crimlab/reference/asciifb.md)
  and
  [`wrtxtfb()`](https://rootcoder007.github.io/crimlab/reference/wrtxtfb.md)
  transliterate text to plain ASCII by a fixed table, keep valid UTF-8
  as it is, and write text as UTF-8 bytes whatever the locale. They are
  in `clab_jrun()`, the command line and Python.

- The summary statistics of the compiled core are exported: the mean,
  variance, correlation and standard deviation
  ([`cmn()`](https://rootcoder007.github.io/crimlab/reference/cstat.md),
  [`cvari()`](https://rootcoder007.github.io/crimlab/reference/cstat.md),
  [`ccorr()`](https://rootcoder007.github.io/crimlab/reference/cstat.md),
  [`csdev()`](https://rootcoder007.github.io/crimlab/reference/cspread.md),
  [`cdistm()`](https://rootcoder007.github.io/crimlab/reference/cspread.md)),
  the normal density
  ([`cnormpdf()`](https://rootcoder007.github.io/crimlab/reference/cnorm.md),
  [`cnormlpdf()`](https://rootcoder007.github.io/crimlab/reference/cnorm.md)),
  moments
  ([`cmom()`](https://rootcoder007.github.io/crimlab/reference/cmom.md)),
  quantiles and robust spread
  ([`cqnt()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
  [`cmed()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
  [`cmad()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
  [`ciqr()`](https://rootcoder007.github.io/crimlab/reference/crob.md),
  [`ctkyfnc()`](https://rootcoder007.github.io/crimlab/reference/crob.md)),
  trimmed and winsorized means
  ([`ctrmmean()`](https://rootcoder007.github.io/crimlab/reference/ctrm.md),
  [`cwinsmean()`](https://rootcoder007.github.io/crimlab/reference/ctrm.md)),
  weighted moments
  ([`cwtd()`](https://rootcoder007.github.io/crimlab/reference/cwtd.md)),
  ranks
  ([`cspr()`](https://rootcoder007.github.io/crimlab/reference/crank.md),
  [`cmidrk()`](https://rootcoder007.github.io/crimlab/reference/crank.md)),
  the covariance matrix
  ([`ccovm()`](https://rootcoder007.github.io/crimlab/reference/ccovm.md)),
  bootstrap means
  ([`cbootmean()`](https://rootcoder007.github.io/crimlab/reference/cbootmean.md)),
  trimmed inverse- probability weights
  ([`cipw()`](https://rootcoder007.github.io/crimlab/reference/cipw.md)),
  the incomplete gamma function
  ([`cgamcdf()`](https://rootcoder007.github.io/crimlab/reference/cgamcdf.md))
  and splitmix64 uniforms
  ([`cunif()`](https://rootcoder007.github.io/crimlab/reference/cunif.md)).
  They are in `clab_jrun()`, the command line and Python.

- [`prechk()`](https://rootcoder007.github.io/crimlab/reference/prechk.md)
  checks a table against a preset before
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  reads it, and names in one error every column the preset reads
  (directly, in a filter or in a recipe) that the table lacks, every
  label outside the codes the publisher lists for its column, and every
  label a mapping recode without a default does not cover. The codes of
  the Ontario segregation, restrictive confinement and deaths-in-custody
  tables are `inst/presets/codes.tsv`, from the ministry’s data
  dictionary; a preset of one’s own takes its codes as a data frame. The
  presets describe the layouts of published tables, and crimlab ships
  none of the tables.

- The OTIS treatment-outcome presets of
  [`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md)
  (`otis-a01-cbatt`, `otis-a01n-cbatt`, `otis-pa-cbatt`,
  `otis-pb-cbatt`, `otis-pc-cbatt`) are tested against rmorie’s pair
  builders on a simulated placement table, with the values stored in
  `tests/testthat/ref/preset.txt` by `dev/ref/preset.R`.

- [`hwpwr()`](https://rootcoder007.github.io/crimlab/reference/hwpwr.md)
  reads the power and energy of the machine’s devices over an interval
  from the operating system’s own counters: RAPL on Linux, the NVIDIA
  management library, the amdgpu and Intel GPU hwmon files, the IOReport
  energy counters of Apple silicon, the Energy Meter Interface on
  Windows, and the battery while it discharges. Each figure says how it
  was obtained; an estimate carries its range, and a device that cannot
  be measured is listed with the reason.

- [`mtrbnd()`](https://rootcoder007.github.io/crimlab/reference/mtrbnd.md)
  bounds the mean outcome under each of ordered sentences and the effect
  of each against the mildest, under the monotone treatment response and
  selection assumptions of Manski and Pepper. It is in `clab_jrun()`,
  the command line and Python.

- [`hsacc()`](https://rootcoder007.github.io/crimlab/reference/hsacc.md)
  scores a hot-spot map on later crimes: the hit rate, the predictive
  accuracy index, the predictive efficiency index and the recapture rate
  index. It is in `clab_jrun()`, the command line and Python.

- [`apcie()`](https://rootcoder007.github.io/crimlab/reference/apcie.md)
  fits the age-period-cohort model of a table of counts by age and
  period with the intrinsic estimator of Yang, Fu and Land: the Poisson
  log-linear fit projected off the null vector of the design. It is in
  `clab_jrun()`, the command line and Python.

- [`spcure()`](https://rootcoder007.github.io/crimlab/reference/spcure.md)
  fits the split-population (cure) model of Schmidt and Witte: a
  logistic share who would ever fail and a Weibull, exponential or
  lognormal time to failure among them, with covariates in both parts.
  It is in `clab_jrun()`, the command line and Python.

- [`lnkpr()`](https://rootcoder007.github.io/crimlab/reference/lnkpr.md)
  links crimes into series: the distance, time and Jaccard similarity of
  behaviour of every pair, a logistic regression of being linked on them
  over the solved pairs, the area under the ROC curve of each, and the
  chance of a link for every pair. It is in `clab_jrun()`, the command
  line and Python.

- [`prohs()`](https://rootcoder007.github.io/crimlab/reference/prohs.md)
  forecasts the risk of each grid cell from recent nearby incidents,
  weighted down with distance and age (Bowers, Johnson and Pease), and
  scores its hot cells on later incidents. It is in `clab_jrun()`, the
  command line and Python.

- [`wdqx()`](https://rootcoder007.github.io/crimlab/reference/wdqx.md)
  gives the weighted displacement quotient of Bowers and Johnson, with
  the counts expected after and the crimes prevented in the treated area
  and its buffer. It is in `clab_jrun()`, the command line and Python.

- [`rptvic()`](https://rootcoder007.github.io/crimlab/reference/rptvic.md)
  measures repeat victimisation: incidence, prevalence and
  concentration, the victims by number of victimisations against a
  Poisson count over the targets at risk, and the time between a
  target’s successive incidents. It is in `clab_jrun()`, the command
  line and Python.

- [`lnroc()`](https://rootcoder007.github.io/crimlab/reference/lnroc.md)
  gives the ROC curve of eyewitness identifications by confidence, the
  partial area with a bootstrap comparison of procedures, the
  diagnosticity ratio and the confidence-accuracy characteristic;
  [`lnsize()`](https://rootcoder007.github.io/crimlab/reference/lnsize.md)
  the fairness of a lineup from mock witnesses (functional size,
  Malpass’s and Tredoux’s effective sizes, the suspect’s share). Both
  are in `clab_jrun()`, the command line and Python.

- [`ibench()`](https://rootcoder007.github.io/crimlab/reference/ibench.md)
  benchmarks each officer’s stops against other officers’ stops in the
  same circumstances (Ridgeway and MacDonald), by propensity weights,
  with the Benjamini-Hochberg false discovery rate over the officers. It
  is in `clab_jrun()`, the command line and Python.

- [`thrtst()`](https://rootcoder007.github.io/crimlab/reference/thrtst.md)
  tests discrimination in searches three ways: the hit-rate test of
  Knowles, Persico and Todd, the threshold test with the discriminant
  distributions of Pierson, Corbett-Davies and Goel fitted by maximum
  likelihood, and a pairwise form of the rank-order test of Anwar and
  Fang. It is in `clab_jrun()`, the command line and Python.

- [`gmixlr()`](https://rootcoder007.github.io/crimlab/reference/gmixlr.md)
  now follows EuroForMix’s model where they differed: only the alleles
  seen are kept apart, every other allele (those of the reference
  profiles not seen among them) is one allele that does not stutter,
  sized as the largest size listed at the marker, and the degradation
  slope may exceed 1. Its log likelihood equals EuroForMix’s to 1e-6 at
  given parameters on eight cases and at the maximum on seven (one and
  two contributors, drop-out, drop-in, stutter, degradation, coancestry;
  `dev/ref/efm.R`); the rest is EuroForMix’s single precision.

- [`cmpsc()`](https://rootcoder007.github.io/crimlab/reference/cmpsc.md)
  now reproduces cmpsR exactly: a segment widened at the ends of the
  signature is cut, not moved; a lag counts when enough points overlap
  (so a segment may hang over the end of the reference); every confirmed
  peak is a candidate; ties are broken as cmpsR breaks them; and lags
  may fall on half points. `full = TRUE` keeps the earlier rule of lags
  wholly inside the reference. The scores, segments, lags and refusals
  equal cmpsR’s on 16 stored cases, edge cases among them.

- [`pedlr()`](https://rootcoder007.github.io/crimlab/reference/pedlr.md):
  the likelihood ratio of two pedigrees, given as tables of people and
  their parents, from DNA profiles at autosomal markers. The people not
  typed are summed over exactly by variable elimination on the alleles
  each received from each parent, so loops, inbreeding and any number of
  generations are allowed; mutation follows the equal or the stepwise
  model, with one rate or separate maternal and paternal rates. It
  equals
  [`kinidx()`](https://rootcoder007.github.io/crimlab/reference/kinidx.md)
  and
  [`patidx()`](https://rootcoder007.github.io/crimlab/reference/patidx.md)
  where they apply, the enumeration of every allele in its tests, and
  pedprobr’s likelihoods.

- [`gmixlr()`](https://rootcoder007.github.io/crimlab/reference/gmixlr.md):
  the likelihood ratio of two accounts of a DNA mixture by the gamma
  model of peak heights, with drop-out, drop-in, backward stutter,
  degradation and coancestry, the unknown contributors summed over and
  the parameters fitted by maximum likelihood under each account, or
  fixed with `par`. Its likelihood is checked against enumeration in R;
  the help page states the model’s assumptions and limits.

- [`cmpsc()`](https://rootcoder007.github.io/crimlab/reference/cmpsc.md):
  the congruent matching profile segments score of two bullet
  signatures, with the multi-segment-length peaks of Chen et al. (2019).
  It follows the algorithm written out in R in its tests and is within
  one segment of cmpsR on simulated signatures.

- [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  reaches the command line, the Python package and the bindings through
  `clab_jrun()`, which does its recodes, filters and recipes in C++ with
  R’s messages and returns the generic’s arguments as `args` with the
  attributes of the R result beside them. Its preset table is
  `cli/pre.inc`, made from `inst/presets/presets.tsv` by `cli/pre.sh`;
  CI fails when the two differ. `cli/test.tsv` holds a case for every
  preset. Outside R, text is ordered by its bytes (R’s order under
  `LC_COLLATE=C`) and lowered as ASCII.

- `clab_jrun()`, and with it the command-line program, the Python
  package and the bindings, now reaches every exported function: counts,
  rates and tables, the audits of risk scores, judges’ leniency, long
  stays, published figures and bands, stocks and flows, the Hawkes fits
  and other events in time, difference in differences, synthetic
  controls and the event study, surveys, few-cluster inference, the
  double machine learner and the super learner, and the spatial methods
  of near repeats, scans, densities, Moran’s I and spatial regression. A
  warning R gives is returned in the result as `warnings`, an array of
  its messages.
  [`lrpost()`](https://rootcoder007.github.io/crimlab/reference/lrpost.md),
  [`impfpr()`](https://rootcoder007.github.io/crimlab/reference/impfpr.md)
  and
  [`lrcal()`](https://rootcoder007.github.io/crimlab/reference/lrcal.md)
  there now give one value as a number, not an array of one, as R does.
  `src/dst.cpp` adds the chi-square (central and noncentral), t, F,
  gamma, beta, Poisson and binomial distribution functions, upper tails
  and quantiles the front ends need; the R functions keep calling R’s
  own. The test of Ctrl-C in a running session measures the stop from
  the signal and sets the signal only once the call’s handler is in
  place, so a loaded machine no longer fails it.

- The command-line program, the Python package and the Go, Rust, Java
  and JavaScript bindings now reach most of the package’s functions
  through `clab_jrun()`: fairness audits, trace evidence, aoristic
  weights and digit tests, sensitivity bounds, trajectory models,
  networks, the reporting checks, the spatial and hot-spot methods,
  spatial weights, toxicology, deterrence and severity, concentration
  and the dark figure, crosswalks and taphonomy, the learners, risk
  terrain and the life course, covariate balance, judge designs and
  Hawkes processes. Each reads its arguments with the R function’s
  checks and words and returns R’s result as JSON; `cli/test.tsv` holds
  the values of the R package for them (`dev/clref.R`). The normal
  quantile and distribution function the front ends need are kernels of
  their own (`src/dst.cpp`).

- The tests no longer call other packages: the values of circular,
  DescTools, geosphere, geepack, igraph, MASS, sandwich, sf, spatstat,
  spdep, rmorie and rmoriebricklayer they compare with are stored in
  `tests/testthat/ref`, with the versions that computed them, by the
  scripts in `dev/ref`. Suggests now holds only testthat, knitr and
  rmarkdown. The `ref` workflow runs those scripts again every week on
  Linux (x86-64 and arm64), macOS (arm64 and Intel) and Windows, and
  fails when a value its platform’s packages give differs from the
  stored one by more than the tolerance of the test that uses it.

- [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  stacks the learners of the package (least squares or logistic
  regression, the lasso, ridge and elastic net, the random forest) into
  a super learner: K-fold cross-fitted predictions, weights by
  non-negative least squares that sum to one, or the discrete super
  learner (the learner of least cross-validated risk), for numeric and
  0/1 outcomes, the folds drawn from a seed and returned. With one
  learner the stack is that learner.
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
  takes it as a nuisance learner (`lrn = "stack"`, with `stk`, `nstk`
  and `disc`).

- Ctrl-C (Esc in the R GUIs) now stops every long computation: the
  resampling, permutation, bootstrap and simulation loops, the iterative
  fits (elastic net, logistic and count regressions, EM, quasi-Newton,
  Weber point, k-means, power iterations, Louvain, synthetic control
  weights, fixed-effects sweeps) and the scans over pairs check for an
  interrupt every few milliseconds of work, free what they hold and end
  with R’s interrupt, never with a half-made result. In ‘Python’ the
  same loops raise `KeyboardInterrupt`.

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
  sandwich of `sandwich::vcovCL()`), and from an exchangeable GEE (as
  ‘geepack’ estimates it). `proof("cntreg")` names the Lean theorems on
  swapping the groups.

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

- [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  refuses a CSV with a quote left open in its first rows as a file that
  “could not be read cleanly”; it used to call the file empty. A file
  without a newline at its end no longer lets a warning escape from the
  reading of its header.

- [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  ends its refusal of negative exposure values with the word it was
  missing: “a concentration cannot be negative”.

- [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  refuses a gzip file cut short inside its data as “truncated or
  damaged”; it used to read the part before the cut without a word, so
  the report rested on some of the rows only. The data read are now
  checked against the length and CRC-32 at the end of the file.

- [`missgrid()`](https://rootcoder007.github.io/crimlab/reference/missgrid.md)
  keeps each letter of the header under its own column when the names
  have different lengths; the letters of a longer name used to slide
  left under the column of a shorter one. The command line had the same
  fault and is fixed with it.

- [`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)
  refuses an AR(1) fit whose coefficient reaches 1 or beyond with “the
  fitted ARMA part is not stationary.”, as it does for higher orders; it
  used to go on with an infinite or negative starting variance.

- A number with no bounds given badly is refused as “must be one
  number.” (`beta` of
  [`nomutl()`](https://rootcoder007.github.io/crimlab/reference/nomutl.md),
  [`nompr()`](https://rootcoder007.github.io/crimlab/reference/nompr.md)
  and
  [`nomll()`](https://rootcoder007.github.io/crimlab/reference/nomll.md);
  `c0`, `d0`, `e0`, `a0` and `b0` of
  [`dynirt()`](https://rootcoder007.github.io/crimlab/reference/dynirt.md))
  or “must be one number, neither infinite nor missing.” (`rho` of
  [`dtrid()`](https://rootcoder007.github.io/crimlab/reference/dtrid.md),
  `a` of
  [`inspg()`](https://rootcoder007.github.io/crimlab/reference/inspg.md),
  `cap` and `floor` of
  [`bandx()`](https://rootcoder007.github.io/crimlab/reference/bandx.md));
  the messages used to name a bound of “-inf” or “-Inf” that does not
  exist. The command line gives the same messages.

- [`pedlr()`](https://rootcoder007.github.io/crimlab/reference/pedlr.md)
  refuses a pedigree whose tables would pass its size bound before it
  builds them: an untyped child of untyped parents at a locus of 2000
  alleles used to try to allocate 8e9 cells and abort R. The mutation
  matrices are made only with a mutation model, and refused past the
  same bound (more than 8192 alleles); they were made even with
  `mut = "none"`.
  [`gmixlr()`](https://rootcoder007.github.io/crimlab/reference/gmixlr.md)
  no longer lists the genotypes of unknown contributors when there are
  none, and
  [`ibench()`](https://rootcoder007.github.io/crimlab/reference/ibench.md)
  refuses more than 8190 columns of `x`, whose work arrays grow with the
  square of their number.

- The Bayesian samplers
  [`bamscl()`](https://rootcoder007.github.io/crimlab/reference/bamscl.md),
  [`bmdsf()`](https://rootcoder007.github.io/crimlab/reference/bmdsf.md),
  [`bunfld()`](https://rootcoder007.github.io/crimlab/reference/bunfld.md),
  [`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md),
  [`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md)
  and
  [`ordirt()`](https://rootcoder007.github.io/crimlab/reference/ordirt.md)
  refuse an `nsamp` whose draws would pass 2^26 values before any is
  kept, and all seven, with
  [`dynirt()`](https://rootcoder007.github.io/crimlab/reference/dynirt.md),
  refuse a `burn` plus `nsamp` (times `thin`) past the largest integer;
  a large `nsamp` used to abort R, and the sweep count overflowed.
  [`optcls()`](https://rootcoder007.github.io/crimlab/reference/optcls.md),
  [`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md),
  [`nmdsf()`](https://rootcoder007.github.io/crimlab/reference/nmdsf.md),
  [`mlsmu()`](https://rootcoder007.github.io/crimlab/reference/mlsmu.md),
  [`smunf()`](https://rootcoder007.github.io/crimlab/reference/smunf.md),
  [`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md),
  [`indscl()`](https://rootcoder007.github.io/crimlab/reference/indscl.md),
  [`dwnom()`](https://rootcoder007.github.io/crimlab/reference/dwnom.md)
  and
  [`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md)
  refuse an `ndim` whose configurations would pass 2^26 values,
  [`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md)
  more dimensions than legislators,
  [`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md)
  too many replicates and
  [`anchv()`](https://rootcoder007.github.io/crimlab/reference/anchv.md)
  too many categories, all before allocating;
  [`nombt()`](https://rootcoder007.github.io/crimlab/reference/nombt.md)
  no longer keeps its replicates (it draws them twice from the seed, for
  the means and then the deviations), so any `nboot` is only long;
  [`smacf()`](https://rootcoder007.github.io/crimlab/reference/smacf.md)
  can now be interrupted. The command line gives the same messages.

- [`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md)
  refuses a horizon `h` past 2^26 and
  [`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md)
  a grid past 2^26 points before making them (both used to abort R); the
  forecast and the steps of
  [`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md)
  can now be interrupted. The command line gives the same messages.

- [`bkrig()`](https://rootcoder007.github.io/crimlab/reference/bkrig.md)
  counts the points of its blocks before it lays them out and refuses
  more than 2^26: `nq = .Machine$integer.max` overflowed the count and
  crashed R. The block integrals can now be interrupted within a block.
  The command line gives the same message.

- [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  refuses `leads` plus `lags` past the number of rows, which cannot be
  estimated (each event time needs a row), before it builds the
  indicators; a large window used to abort R.
  [`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
  with an `mtch` larger than the number of comparison units matches as
  that many does, instead of sizing its tables by `mtch`. The command
  line does the same.

- [`hwksim()`](https://rootcoder007.github.io/crimlab/reference/hwksim.md)
  refuses an `nmax` past 2^26 and
  [`sdbsim()`](https://rootcoder007.github.io/crimlab/reference/sdbsim.md)
  more than 2^26 `steps`, before the arrays they size are made (both
  used to abort R).
  [`kscan()`](https://rootcoder007.github.io/crimlab/reference/kscan.md)
  compares each simulated maximum as it is drawn instead of keeping all
  `nsim` of them, so a large `nsim` is only long, and can be
  interrupted. The command line gives the same messages.

- [`lnroc()`](https://rootcoder007.github.io/crimlab/reference/lnroc.md)
  no longer keeps its bootstrap replicates: it draws them twice from the
  seed, summing the areas and then their squared deviations, so the
  errors are those of the kept replicates and a large `nboot` (which
  used to abort R) is only long.

- [`inlhist()`](https://rootcoder007.github.io/crimlab/reference/inlhist.md)
  refuses more than 2^26 `bins`; the command line gives the same
  message.

- [`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md)
  gives the negative binomial bound for large counts (more than about
  750 a period in the reference window): the search started where the
  probability underflows to 0 and never ended. Counts so large that no
  bound can be found are refused.
  [`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md)
  stops selecting when no eigenvector has a finite criterion (values
  whose sums overflow); it used to add one that does not exist and abort
  R.

- [`dualcont()`](https://rootcoder007.github.io/crimlab/reference/dualcont.md),
  [`frmaud()`](https://rootcoder007.github.io/crimlab/reference/frmaud.md)
  and the readers of
  [`cnsprof()`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md)
  results check the columns they read even when the class says where the
  object came from, once their other arguments are checked; an object of
  that class with a column missing crashed R.

- [`tpdmc()`](https://rootcoder007.github.io/crimlab/reference/tpdmc.md)
  refuses more than 2^26 `steps` before it makes the path; the command
  line gives the same message.

- [`sirex()`](https://rootcoder007.github.io/crimlab/reference/sirex.md)
  checks `level` before it reads it; an empty `level` was read past its
  end and crashed R.

- [`tfrp()`](https://rootcoder007.github.io/crimlab/reference/tfrp.md)
  refuses a `ymax` of 2^26 or more and
  [`aorwt()`](https://rootcoder007.github.io/crimlab/reference/aorwt.md)
  more than 2^26 bins, before anything is made;
  [`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md)
  makes as many records as asked, or stops with “fairsim: not enough
  memory for `n` values.” The command line gives the same messages.

- [`crwalk()`](https://rootcoder007.github.io/crimlab/reference/crwalk.md)
  makes its simulated paths once, as the matrices it returns (they were
  made twice), and stops with “not enough memory for `h` times `nsim`
  values.” when they cannot be made;
  [`vgemp()`](https://rootcoder007.github.io/crimlab/reference/vgemp.md)
  refuses more than 2^26 distance bins. The command line does the same.

- [`lnsize()`](https://rootcoder007.github.io/crimlab/reference/lnsize.md)
  refuses a lineup size `k` past 2^26; the command line gives the same
  message.

- [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  refuses more than 2^26 `draws` for the sensitivity interval, and
  [`wcbx()`](https://rootcoder007.github.io/crimlab/reference/wcbx.md) a
  `reps` whose draws (one value per cluster and outcome) would pass
  2^26, in R and in the command line. The draws of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md),
  the bootstrap of
  [`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md),
  the placebo replicates of
  [`sdidx()`](https://rootcoder007.github.io/crimlab/reference/sdidx.md)
  and the permutation and power replicates of
  [`wcbx()`](https://rootcoder007.github.io/crimlab/reference/wcbx.md)
  can now be interrupted.

- [`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md),
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
  and
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  refuse an `nlam` past 2^26, in R and in the command line.

- [`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md)
  refuses a `width` that divides 360 into more than 2^26 sectors, or
  into more columns (areas times sectors) than an integer counts; a tiny
  `width` used to overflow the count of columns. The command line gives
  the same messages.

- The robustness iterations of
  [`lwss()`](https://rootcoder007.github.io/crimlab/reference/lwss.md),
  the replicator steps of
  [`insphs()`](https://rootcoder007.github.io/crimlab/reference/insphs.md),
  the cells of
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md),
  the EM steps of
  [`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md)
  and the bootstrap refits of
  [`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md)
  can now be interrupted. Loops and sizes that multiplied two integer
  counts are counted in a wide type.

- [`exprun()`](https://rootcoder007.github.io/crimlab/reference/exprun.md)
  and
  [`runfx()`](https://rootcoder007.github.io/crimlab/reference/exprun.md)
  refuse a `cap` past 2^26,
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md)
  an `ngrid` past 8192 (its grid has `ngrid` squared points) and
  [`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md)
  an `nq` past 2^26. These used to stop with R’s own “cannot allocate
  vector” error. The command line gives the same messages where it has
  the function. The cluster bootstrap of
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  and
  [`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)
  makes its draws a block at a time instead of holding the clusters by
  `nboot` counts at once, with the same draws in the same order, so any
  number of clusters and draws fits.

- [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
  and
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  refuse more inner folds (`nfin`) than rows of `x`, before their work
  arrays are sized by it; the command line gives the same message.

- [`cconc()`](https://rootcoder007.github.io/crimlab/reference/cconc.md)
  refuses an `nsim` past 2^26,
  [`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md)
  an `nboot` past 2^25 and
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md)
  more than 2^26 colours, instead of R’s “cannot allocate vector” error;
  [`cbootmean()`](https://rootcoder007.github.io/crimlab/reference/cbootmean.md)
  and
  [`cunif()`](https://rootcoder007.github.io/crimlab/reference/cunif.md)
  make results of any length R can hold, and stop with “not enough
  memory for `B` values.” (or `n`) when it cannot. The command line
  gives the same messages.

- [`stacel()`](https://rootcoder007.github.io/crimlab/reference/stacel.md)
  with a `top` above the cells of its grid keeps them all, as that many
  does, instead of sizing its work arrays by `top`.

- No exception leaves a C entry point of the library: every `clab_*`
  function returns `CLAB_ENOMEM` (new, 4) when an allocation fails in
  it, `clab_jrun()` then answers “fn: not enough memory for this size.”,
  and the command-line program prints that and exits with status 1 for a
  failed allocation anywhere, never aborting. Python raises a
  `MemoryError` with the message. The WebAssembly build, which has no
  exceptions, stops the module and the JavaScript binding reports the
  same words (`clab_oomed()`, new) and asks for the module to be loaded
  again. In R the entry points still give the function’s name.

- Sizes and products of sizes are 64-bit everywhere they were `long`,
  which is 32 bits on Windows: the work arrays of the learners
  ([`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md),
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md),
  [`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md)),
  the neighbour tables, the bins of
  [`rddens()`](https://rootcoder007.github.io/crimlab/reference/rddens.md)
  and others. `clab_dmlxsz()` and `clab_stksz()` now return `int64_t`
  and take an `int64_t *`.

- The quantile functions of the command line, the JSON entry and Python
  give R’s answers at a point mass: the beta quantile with both shapes 0
  is 1/2 at p = 1/2 (it was 1) and compares p itself in the upper tail
  as R’s [`qbeta()`](https://rdrr.io/r/stats/Beta.html) does; the gamma
  and chi-square quantiles with shape or degrees of freedom 0 are
  infinite at p = 1 (they were 0).

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
  corrections of `spatstat.explore::density.ppp()`) with the smallest
  areas holding each share of the incidents;
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
  (`spdep::moran.test()`, `spdep::moran_bv()`), one variable or a sweep
  of columns;
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
