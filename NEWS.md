# crimlab 0.0.4

* Fixed: `spill` from the command-line program, and from the Python, Go,
  Java, JavaScript and Rust packages, read outside a list when the fit it
  was given named more lagged terms than it had terms (a fit not made by
  `sarml`, as text). Such a fit is now refused. The R function was not
  affected. Found by the fuzzing of every change.
* No function changes. The command-line program and the library build on
  systems with an older C library (glibc before 2.34, where the loader and
  the threads are outside libc): the Makefile and the Python, Go and Rust
  builds name those libraries at the link. A build with an older compiler
  (gcc 8) no longer stops at a warning, and `make install` creates each
  directory in a command of its own.
* A Windows arm64 archive of the command-line program, beside the x86_64
  one; `install.ps1` and the Scoop manifest take the one that fits the
  machine.
* The pages of the package on PyPI, npm, crates.io, pkg.go.dev and Maven
  Central carry the badges of the main README, and the README lists every
  way to install.
* The source archive no longer carries build output of the Python package,
  which the archive of 0.0.3 included by mistake.
* The Linux wheels on PyPI are a tenth of their size: the extension was
  built with debug records (about 40 MB a wheel at 0.0.3), and is now built
  without them. A wheel over 20 MB stops the build.
* A notice of a newer release. In an interactive R session, attaching the
  package says in one line when a newer release exists, and `import
  crimlab` does the same at the Python prompt or in a notebook: the page
  of the latest release is asked for at most once a day, with a limit of
  two seconds, and a script or a test never asks.
  `options(crimlab.update = FALSE)` in R, or the environment variable
  `CRIMLAB_NO_UPDATE_CHECK`, turns it off. The command-line program asks
  only when told to, with `crimlab --update`, and keeps nothing. Nothing
  is sent but that request.

# crimlab 0.0.3

* No function changes. The command-line build can lower the compiler
  flags of its JSON layer alone (`JOPT`), which the Fedora package does
  on riscv64: that one file took the emulated builders over an hour, and
  the packages of 0.0.2 were not built there.

# crimlab 0.0.2

* The fits of `ppmlfe()`, `ivfe()`, `ppmlcf()`, `expgmm()`, `fepnl()`,
  `gmmpnl()`, `didfe()`, `rdloc()` and `sarml()` answer to the usual
  accessors: `coef()`, `confint()`, `vcov()` and `nobs()` (help page
  `racc`). The intervals are the fit's own, at its own level.

* A file written under a name with no suffix is given one: `wrtxtfb()`
  adds `.txt`, and `yoycsv()`, `yoytsv()`, `yoyjs()`, `yoymd()`,
  `yoyhtml()` and `yoypdf()` the suffix of their format. Each returns
  the path it wrote.

* The rOpenSci standards notes describe the package as it is now: the
  regression standards are documented (`src/srrre.h`), the general and
  spatial notes written when the package was small are brought up to
  date, every internal R function has a roxygen comment, and
  `tests/README.md` says how the full test suite is switched on.

* The argument that names a method is `mthd` (it was `meth`), and so is
  the part of a result that records it: in `lodimp()`, `stgdid()`,
  `wcbx()`, `sarml()`, `cortbl()`, `topcors()`, `evfilt()`, `scboot()`,
  `chomap()`, `smmap()` and `co2fp()`, in R and in every other language.
  A call that names the old argument is refused.

* Plots are black and grey on a white ground with a black frame and a
  bold title, as a journal prints them: the tinted panel and its grid are
  no longer the default. `plotopt(palette = "okabe")` brings the colours
  back, and `plotopt(grid = TRUE, panel = "#EEF1F6")` the panel and grid.
  The `"grey"` palette starts at black and has no grey too pale to read
  on white.

* Intervals are drawn with a cap at each end. `plot()` of `evtst()` and
  `stgdid()` with `pal = "diverging"` is blue before adoption and red
  after. `plot()` draws the power curve of `wcbx(mthd = "power")`.
  `plotevgrid(x, style = "stata")`, and `plot()` of `evtst()` and
  `stgdid()` with it, draws hollow points joined by a line with dashed
  intervals. The Goodman-Bacon panel of `plot(didfe())` names the two-way
  estimate its line marks.

* New functions for the review of restrictive confinement, after the
  analyses of Sprott and Doob (2021) and of Sprott, Doob and Iftene (2021)
  of Canada's Structured Intervention Units, in R and in every other
  language. `rvwtim()`: whether a review due by a set day was held, the
  stays past that day never reviewed, and how long the reviews held took.
  `entcmp()`: the days of a stay on which an entitlement (hours out of
  cell, human contact) was missed, by how often, with the hours received
  on those days and the part refusals account for. `dmvar()`: whether
  decision makers decide alike, with each one's rate, the chi-square test,
  and the test with the small ones and with each one left out.

* New plots. `plot()` of `wcbx(mthd = "perm")` draws the permuted
  statistics of each outcome with the observed one, as a histogram or a
  density (the result now keeps them, as `null`, in every language).
  `plot(didfe(loo = TRUE), type = "hist")` draws the leave-one-out
  estimates as a histogram with a summary. `plotcoh()` draws the mean
  outcome of each adoption cohort around its adoption, and returns the
  points in every language. `plot()` of
  `yoypx()` with `what = "value"` draws the values over the periods.

* The command-line program and the library are built with the arithmetic
  of SSE2 on 32-bit x86, as every other platform computes. With x87
  arithmetic and a distribution's options (link-time optimisation, no
  frame pointer) three results of the test table differed in their last
  digits, and the openSUSE i586 package of 0.0.1 failed its checks. A
  build that names `-mfpmath` itself keeps its choice.

* A check of CRAN runs one test file, the compiled self-test of every
  kernel and one call of every exported function; the whole suite, which
  takes minutes, runs in the project's own checks on every push
  (`NOT_CRAN=true`).

* The package installs faster: R compiles its C++ sources in 40 groups
  (a unity build, `src/u_*.cpp`) and no longer one file at a time, so the
  headers are read 40 times, not 337. Nothing a caller sees changes.

* Outside R the plots are drawn too. `plotsvg()` (Python, Go, Rust,
  Java, JavaScript and the command line) runs a plot function on its
  arguments and returns the picture as a standalone SVG document with the
  rows drawn: `plotlrnz()`, `plotgain()`, `plotcdfs()`, `plotrdg()`,
  `plotvln()`, `plothex()`, `plotfnl()`, `plotcal()`, `plotmult()`,
  `plotbump()`, `plothrz()`, `plotevgrid()` and `plotcoh()`, and every
  result R draws with `plot()` (`evtst()`, `stgdid()`, `synthc()`,
  `sdidx()`, `didfe()`, `cntreg()`, `kmrec()`, `rdloc()`, `indgap()` and
  the rest, 109 plots in all; `bandx()` is left to R, since its plot needs
  a statistic given as an R function), with the choices `type`, `what`,
  `units` and `terms` of those plots and the look
  of the R plots (black and grey, a black frame, a bold title; `palette =
  "okabe"` for colours, `grid = TRUE` for a grid). One drawing code in the
  C++ core serves every language.

* `plot()` of `synthc()` and `sdidx()` draws the treated path solid and
  its comparison dashed, as plain lines, with a dotted line where the
  policy was adopted. `plot(x, type = "o")` puts a point mark on every
  period, as before.

* Every package carries the title of the R package, "Statistical Toolkit
  for Criminology and Sociolegal Research": the Python, JavaScript, Rust,
  Java and Go packages, the command's manual page, the C header and the
  Homebrew, Scoop, Debian, RPM and container descriptions said
  "criminology and sociolegal statistics". The README of the Python,
  JavaScript and Rust packages, which is what each registry shows, says
  how to cite crimlab.

* Fixed: `ecolrg(family = "negbin")` on counts with no overdispersion.
  The search for the size sends it off without bound there, and past
  about ten million its score is no more than rounding; on 32-bit x86
  without SSE2 (x87 arithmetic) a step taken from that noise ended at a
  wrong fit, with a lower likelihood. The search now ends at the cap
  once the size passes a million, where the fit is the Poisson fit to
  six digits. Results elsewhere do not change.

* The tests of the command run on x87 as well: its test table compares
  the numbers inside a text result, and two cases whose last digits
  differ there have a bound.

# crimlab 0.0.1

* The results of the Indigenous health and wellbeing functions have
  `plot()` methods: gaps between groups with their intervals
  (`indgap()`, `trchlthgap()`), wellness beside distress (`dualcont()`,
  `welprof()`, `holbal()`), community indices and Community Well-Being
  scores (`commidx()`, `cnscwb()`, `cwbscr()`), standardised prevalence
  (`stdprev()`), ordered ratings (`ordcmp()`) and change between surveys
  (`healchg()`). A dot plot keeps every row label, in letters sized to
  the rows, and leaves room for its key above the data; lines over
  periods written as text ("2020/21") are drawn in the order given and
  labelled as written.

* Open data: `odtab()` loads a table of the project's open data service by
  its key and `odls()` lists the tables, with their rows, columns and
  SHA-256. Everything about a table is read from the service's manifest
  at the time of the call, so a table its publisher replaces is loaded as
  it is now, with no new release. A file is used only when its SHA-256 is
  the one the manifest lists; tables are kept in
  `tools::R_user_dir("crimlab", "cache")` while they still match. `cols`
  picks columns without regard to case or underscores, and `base` (or
  `options(crimlab.data)`) points at a mirror.

* Placement volatility and alert complexity: `plcvol()` counts the changes
  of placement region in each person-year (within a record and between
  consecutive records) and `alrtcx()` the distinct alert configurations
  its records span, as defined in Ruhela's study of Ontario's restrictive
  confinement records. `otpyr()` reduces a release of those records (a
  data frame, or a key of the open data service read by `odtab()`) to one
  row per person and fiscal year, finding the columns without regard to
  case or underscores and taking the fiscal years from the data. `cxvol()`
  gives, pooled or for each fiscal year, the descriptives, the crude rate
  ratio, the one-to-one propensity matching, the count model with a random
  intercept per regional cluster and the bounds on the effect without
  ignorability; `cxdml()` the interactive double machine learning
  estimates (ATE and ATTE); `cxyoy()` the year-over-year table of
  regional changes per 1,000 placements and per 1,000 person-years with
  exact intervals. A single fiscal year with a year term, or a sample
  with no complex or no other person-year, is refused in words.
  `cntmix()` fits a Poisson or negative binomial count with a normal
  random intercept per cluster by the Laplace approximation, and
  `cxsim()` simulates the bias, spread and coverage of the estimators,
  with or without an unobserved confounder.

* `ctabx()` bounds Cramer's V right at any total: for a chi-square past
  about 1.5e6 (a 2 x 2 table of 1e8 counts) both bounds in R fell below
  the estimate, with warnings from R's `pchisq()`, and on the command line
  and in Python the bounds took 149 s at 1e9 counts. R, the command line
  and Python now take them from one routine of the core, an integral of
  the noncentral chi-square whose work does not grow with the counts
  (about a millisecond a bound, to 1e-14 of references of 34 to 40
  digits). With `exact = TRUE` the odds ratio of a 2 x 2 table and its
  interval are summed about their mode: a table of 1e6 counts took a
  minute.

* `pwranv()` gives power 1 for a large effect (0.9999999486 at `f =
  1000`) and returns at once for any `f`; `f = 2^31` and larger never
  returned.

* `gstr()` takes the same time whatever the counts: the binomial tail of
  its exposure test summed the terms one by one (20 s for counts near
  2e6 with `nsim = 0`, no end for larger ones).

* `opteffort()` gives finite efforts that add to the total for any total:
  from about 1e4 the multiplier underflowed and the first cell's effort
  was `Inf`.

* `locgi()`, `olsdx()` and `sarml()` give the same results for values
  shifted by a constant (and `locgi()` for values scaled): with values
  large next to their spread `locgi()`'s z-values moved by 0.4 (1e8) or
  became `NaN` (1e9), `olsdx()` moved its slope (1e6) or refused the
  column as collinear (1e7) and `sarml()` refused it (1e6).

* `gtwreg()` fits an integer `x` with `cons = FALSE` (it stopped in R's C
  words).

* `plothrz()` refuses a period that is not finite; with one it never
  returned.

* Refusals and warnings start with the function's name, in the
  package's words: about forty functions refused some inputs in base R's
  words (a date where a column name goes, a function where a vector goes,
  text that is not valid UTF-8, a path in a folder that does not exist, a
  count past what can be allocated, a 1 x 1 matrix for a number), and
  the warnings of `stdprev()`, `dosebin()`, `hetfx()` and `jntx()` did
  not name their function. A text argument that is not valid UTF-8 is
  refused in R as on the command line ("`ref` holds text that is not
  valid UTF-8"); in Python with the same words. `plotgain()` no longer
  says that counts whose total overflows are all 0. A date or a time
  given where the scaling and IRT models take a count (`burn`, `nsamp`,
  `ndim`, `maxit`) is refused; it was read as its days or seconds since
  1970 and ran for hours.

* `dktflw()` takes its clearance intervals from the core's beta quantile,
  as the command line does: R's lost its digits, with warnings, once one
  count was 1e13 times the other.

* Proofs: `inst/lean/Volat.lean` (the two measures: bounds, the cases of
  no change and of one record, relabelling, order of the records with a
  counterexample and the sorted case, adding a record, pooling, weighted
  and unweighted means, and two counterexamples on what alert complexity
  counts), `inst/lean/Ident.lean` (identification by adjustment, weighting
  and augmentation, double robustness, exact matching, the log-link rate
  ratio, a negative control, Manski's bounds and the non-identification of
  the counterfactual arm) and `inst/lean/Custody.lean` (Lakner's average
  daily population, length of stay and expansion estimator with the
  manual's worked examples, and the arithmetic of year-over-year changes,
  including what a first-to-last column means once a year is added).
  `proof()` lists them.

* `odtab()` gives the values the file holds: every column is read as
  text and becomes numbers only when each value is a plain decimal number
  held exactly, so identifiers such as `"007"`, `"7"` and `"0007"` stay
  distinct and a column of `"F"` stays text (it became `FALSE`). A line
  with more or fewer fields than the header is refused, naming the line
  (a long line was wrapped into a row of its own), and so is a table
  whose number of rows is not the one the manifest lists. The manifest
  reader takes all JSON: an empty key, a string of any length (one of a
  few megabytes failed in PCRE), and characters past U+FFFF written as
  surrogate pairs.

* Messages quote invalid UTF-8 bytes of the input as `<xx>`, in every
  refusal and warning of the package and of its C++ binding: `dspbm()`,
  `clprep()` and about twenty others carried the bytes themselves, and
  `substr()` or `grepl()` on such a message failed.

* The keys that join several columns (the person-years of `plcvol()` and
  `alrtcx()`, the units and days of `expjn()`, `smkiv()`, `ccross()` and
  `updwn()`) no longer merge labels holding a carriage return: each value
  is escaped before the values are joined, so no label can imitate the
  join.

* `qltyrcd()` reads its flags once: a vector of a million flags took
  minutes, copying a row per value, and is now read in a fraction of a
  second.

* Arguments of the wrong type are refused in the package's words, not
  base R's: `gtwreg()` takes `x = NULL` (a model with the constant only,
  which stopped on its column names) and refuses a function or an
  environment as `x`; `grdbin()` a `col` that is not one name;
  `xfrvrf()` an `imported` that is not a vector; `dagcol()` a `triple`
  that is not labels; `lrcce()` an `x` that names no column; the plots
  a `col` that is not colours (a function was warned about, a data
  frame failed in the graphics).

* The help of `vrfpol()` says what it refuses with an R error: only a
  `poll` that is not one pollutant name, as the command line refuses it;
  it said the function never stops with an R error.

* `bayout()` gives R's `qnbinom()` bound for every reference sum, in R,
  the command line and Python. For sums near 4,800 (a week of about 800)
  the probability of 0 was a subnormal number: the bound came out too low
  (805 for 856, raising alarms wrongly) and for nearby sums the search
  never ended (`bayout(1:810)` did not return). The quantile is now a
  search of at most about 160 evaluations of the distribution function,
  and the Poisson quantile of the command line and Python takes the same
  search.

* `vrfpol()` is in `clab_jrun()`, the command line and Python, with the
  report R gives: the demo data are R's draws, and the exposure file is
  read as R reads plain and gzip UTF-8 text. A bzip2 or xz file and
  UTF-16, UTF-32 or Latin-1 text, which R decodes, are an error report
  outside R.

* New: `tabrpt()`, one call on a published table: `catchk()` on its
  groups, `yoypx()` on the change between periods, `grprr()` on the rates
  and `trndrr()` on each group's trend.

* New: `ovrchk()` flags the cells of a coded table of oversight reports
  that are not of the form of their field (a real date, a count, a flag)
  or of a pattern such as a case number.

* New: `ovrkwd()` marks the reports that mention the words of keyword
  groups, whole words in any case, with a preset for mental health and
  race and any groups of one's own.

* New: `ovrsoc()` counts the subject officials of an oversight report by
  the first of six rules that applies (a roster, numbered tags, a plural,
  a singular mention, a statement of none) and says which rule and words
  gave the count.

* New: `ovrpar()` codes police oversight reports into the fields of a
  schema (`ovrsch()`), by rules that are deterministic and offline, one
  row per report.

* New: `ovrsch()`, the field schemas of oversight reports: each field
  read by a rule (text after a cue, a number, a date in the cue's
  sentence, numbered tags, a flag and others), with a preset for the
  director's reports of Ontario's Special Investigations Unit.

* New: `ovrdte()` reads the first date written out in a text, in English
  or French, as `YYYY-MM-DD`.

* New: `ovrtxt()`, the plain text of a police oversight report's page,
  with each block on a line of its own so that section headings stand
  alone; the first of the tools that code oversight reports.

* New: `synfix()` draws example data in the shape of a published table:
  presets with the columns of Ontario's corrections tables (segregation,
  deaths in custody, use of force), or any table given as the values of
  its columns. The draws are R's after `set.seed()`, in every language.

* New: `pscred()`, the credit for custody before sentence under s. 719(3)
  and (3.1) of the Criminal Code, at up to one and a half days a day, and
  the time left to serve.
* New: `olsdx()`, least squares with the checks Algahtany, Kumar and
  Barclay (2022) report for their weather and crime regressions
  (standardised coefficients, tolerance and variance inflation, the
  Durbin-Watson statistic, leverage, studentised residuals and Cook's
  distance) and the expected count, or a new month's, at given weather.

* Preset `laniyonu-gent-cnt`: the tract layout of Laniyonu (2018), stops
  on the three-level gentrification factor (the `gent` recipe) with the
  population as the offset and the crime and 311 controls, for
  `cntreg()`.

* New: `dagcol()`, the colliders of a directed acyclic graph given by its
  edges, shielded or not (v-structures), and whether a given triple is
  one.

* New: `tpsurf()`, a thin-plate spline surface through values at points
  in the plane (Duchon 1977), smoothed by `lam`, on a grid with `print()`
  and `plot()` methods.

* New: `eofst()`, the empirical orthogonal functions of a matrix of
  periods by places: spatial patterns, their series over time and the
  share of the variance of each.

* New: `prpmap()`, a proportional symbol map: circles whose areas are
  proportional to the values, drawn largest first.

* New: `ptmap()`, a point map in base graphics, with clusters from
  `dbscn()` or given, each in its own colour and shape with its number at
  its centre; `dbscn()` results now have a `plot()` method that draws
  them through it.

* New: `smmap()`, small multiples of a choropleth map, one panel per
  period with the breaks and colours shared by every period.

* New: `chomap()`, a choropleth map in base graphics from polygon
  coordinates, with no spatial package: quantile, equal or given class
  breaks computed in the compiled kernel.

* New: `gocoef()`, the slopes of a `gologt()` fit laid out by threshold,
  with their standard errors and Wald tests.

* New: `gologt()`, the generalised ordered logit (Peterson and Harrell
  1990; Williams 2006): an intercept and slopes for each threshold of an
  ordered outcome, chosen covariates held parallel, with the Wald test
  of each covariate's slopes and the likelihood-ratio test against
  proportional odds.

* New: `ulcate()`, the U-learner of conditional average treatment
  effects (Kunzel et al. 2019) on the cross-fitted learners of `dmlx()`,
  with a final learner of its own and a clipped denominator.

* New: `xpshap()`, the sampling Shapley values of Strumbelj and
  Kononenko (2014) for single predictions, with their standard errors.

* New: `xpale()`, the accumulated local effects of one input (Apley
  and Zhu 2020), with the binning and centring of their 'ALEPlot'.

* New: `xpcp()`, ceteris paribus profiles: the prediction for single
  rows as one input moves over a grid.

* New: `xppdp()`, the partial dependence of a prediction on one input
  (Friedman 2001) over an even or a quantile grid.

* New: `xpimp()`, the permutation importance of a risk tool's inputs
  (Breiman 2001; Fisher, Rudin and Dominici 2019), for an R function or
  the learners of `suplrn()`, with the shuffles drawn in the compiled
  core and returned.

* New: `mipmm()`, multiple imputation of numeric columns by predictive
  mean matching in chained equations, with an analysis pooled by Rubin's
  rules and the degrees of freedom of Barnard and Rubin, after
  `mcartst()`.

* New: `rscchi()`, the Rao-Scott test of independence in a two-way table
  from a stratified, clustered and weighted survey, first and second
  order and as an F, with the generalised design effects, beside
  `svymx()` and `ctabx()`.

* New: `rcttml()`, targeted maximum likelihood for a randomised trial with
  covariate adjustment (Moore and van der Laan 2009), for continuous and
  binary outcomes, beside the g-computation and the difference in means.

* New: `trnsp()`, a trial's average effect transported to a target
  population, or generalised to the population it came from, by the
  outcome model and by inverse odds weighting, with bootstrap errors.

* New: `gensw()`, the sampling score of a trial against its target
  population and the weights that generalise or transport its units,
  with the logit standardised difference, Tipton's generalisability
  index, the effective size and the covariate balance.

* New: `specrv()`, the specification curve of Simonsohn, Simmons and
  Nelson (2020): the treatment's estimate for every outcome and every
  subset of optional covariates, with the joint inference of the median,
  the count significant and the summed t against permuted or
  bootstrapped null curves.

* New: `plntfx()`, the recovery of planted effects: the outcome shuffled,
  an effect of known size added to the treated and the analysis run
  again, reporting bias, spread, coverage and rejection.

* New: `falsctl()`, falsification with negative-control outcomes: the
  treatment's coefficient on outcomes it cannot affect, with classical or
  HC1 errors, Kennedy's permutation p-values and an adjustment for the
  number of controls.

* New: `pwrsim()`, the power of a two-arm trial by simulation, for normal,
  binary or count outcomes, unequal arms and Welch's or the permutation
  test, with its Monte Carlo error.

* New: `pwranv()`, the power of the one-way analysis of variance from the
  noncentral F distribution, or the group size that reaches a given
  power.

* New: `blkprm()`, the permutation test of a treatment effect with the
  assignment permuted within blocks, by enumeration or by draws, for the
  difference in means or the block-weighted difference.

* New: `rndalc()`, randomisation lists by simple, permuted-block or
  stratified block randomisation from a seed, the same list in R and on
  the command line.

* `plotrdg()` and `plotvln()` take their densities from the compiled
  kernel that the command line and the other languages use, so the curve
  is the same everywhere and on every version of R. The kernel agrees
  with `stats::density()` of R 4.3 and later to rounding.

* New: `cemtch()`, coarsened exact matching (Iacus, King and Porro 2012)
  with the weights for the effect on the treated and the L1 imbalance,
  as Sampson and Winter (2018) match lead-exposed children.

* New: `lrcce()`, the CS-ARDL long-run effects and adjustment speed
  from the unit estimates of `ccepnl()`, as Ditzen (2021) sets them out.

* New: `ifesc()`, the generalized synthetic control of Xu (2017) on the
  interactive fixed effects model of Bai (2009), with staggered adoption,
  the number of factors by cross-validation and a bootstrap of units.

* New: `ppmlcf(family = "nb")`, the control function with a negative
  binomial second stage and a cluster bootstrap of both stages, as
  Feigenbaum and Muller (2016) estimate lead's effect on homicide.

* New: `ivfe(se = "dk")`, the standard errors of Driscoll and Kraay
  (1998), robust to correlation across units and over time, as Batkeyev
  and DeRemer (2023) report them.

* New: `idwint()` carries monitor readings to places by inverse-distance
  weighting of the nearest monitors and averages them by zone, as
  Rahman et al. (2023), Wongchanapai et al. (2025) and Wesselbaum (2022)
  assign exposure; `gtwreg()` is the geographically and temporally
  weighted regression of Kuo and Putra (2021), its bandwidth and time
  scale chosen by cross-validation or AICc; `asrule()` gives the
  association rules between two binned variables that Wongchanapai et
  al. (2025) mine with Apriori.
* The command line computes the quantiles and medians of `jntx()`,
  `ccross(ns =)` and `dosebin()` as R's `quantile()` and `median()` do,
  and the restricted cubic spline's cube as a product in both, so that a
  spline surface of `jntx()` is the same fit in R and outside it (it
  differed by about 1e-8).

* Fixed, after a review of the checks: `ppmlfe()` and `ivfe()` with
  `se = "twoway"` and other than two fixed effects, `rddens()` and
  `rdloc()` with fewer than 10 rows, `rdloc()` with a constant treatment,
  a refused row of `polbrdarea()` and an outcome `polbrdn()` lacks for
  PM2.5 are refused in words that name what is wrong; `tailwd()` refuses
  a negative speed and `dwnexp(lonlat = TRUE)` a latitude past 90;
  `windiv()` gives winds that cancel, and calm, no direction (`NA`); a
  whole number past 2147483647 is refused naming that limit; and the bin
  labels of `invsn()`, `dosebin()` and `updwn()` write `Inf` without a
  space, as the command line does.

* Fixed: Python `orrisk()` refused `or_`, although an argument named as a
  Python keyword takes a final underscore elsewhere (`from_`, `lambda_`).

* Fixed: `inlhist()` lost a value when the range of the values overflowed
  (`c(-1e308, 1e308)`); the breaks are now built from halves there.

* `windiv()` with very narrow sectors builds its column names in less than
  half the time: each edge is formatted once.

* `ccross()` takes a day column of ISO date text or of days since
  1970-01-01, as the command line gives dates, and `smkiv()` reads ISO
  date text as dates, in R and on the command line alike; before, R's
  `ccross()` refused numbers and the command line could not give
  `smkiv()` dates.

* Fixed: `irtem()`, `optcls()`, `dwnom()`, `cjrirt()`, `nomll()`,
  `irtlik()` and `nombt()` accepted votes other than 1, 0 and `NA`, and a
  matrix of no votes; they refuse them now in the words of `anomin()`.
  `ordocl()` refuses a matrix of no responses, and a vector at once.

* Fixed: `clncolnms()` could return repeated names (`c("x", "x", "x_2")`
  gave `"x_2"` twice); a number whose name is taken is now skipped. A tab
  inside a name is a word break, as the help says, and no longer deleted.

* Fixed: `mcartst()` refused columns in different units (income next to
  a rate) as collinear; the rank check and the EM step's eigenvalue floor
  are now on the scale of the correlations.

* Fixed: `rdloc(bw = "ik")` and `rddens()` refused a running variable in
  metres or in thousands of kilometres; the bandwidths now scale with the
  running variable, and the estimates do not change.

* Fixed: `ppmlfe()`, `ivfe()` and the other least-squares fits refused
  columns about 1e7 apart in scale as collinear; a column is now judged
  against its own size, and a column the fixed effects absorb is found by
  an explicit check in `didfe()`, `wcbx()` and the tests of `ivfe()`.

* Fixed: `plot.lnkpr()` returned an ROC curve of `NA` and an area of `NA`
  when some crimes were unsolved; the pairs with an unsolved crime are now
  left out, as `lnkpr()` leaves them out of its fit, and the area is the
  model's AUC.
* `polrpttxt()` is in `clab_jrun()`, the command line and Python: it
  writes the text of a `vrfpol()` report given as JSON, its numbers as
  R's `sprintf()` and `format()` write them. `vrfpol()` itself runs in R
  only.

* `yoypdf()` is in `clab_jrun()`, the command line and Python, writing
  the PDF R writes; without `file` the PDF is returned as text, and
  `yoywr()` there writes PDF too.

* `yoyhtml()` is in `clab_jrun()`, the command line and Python, writing
  the page R writes; without `file` the page is returned as text, and
  `yoywr()` there writes HTML too.

* `prechk()` is in `clab_jrun()`, the command line and Python, with the
  publishers' codes of `codes.tsv` built into the command line as the
  presets are (`cli/codes.sh` writes `cli/codes.inc`, and CI checks it
  is up to date).

* `plotopt()` is in `clab_jrun()`, the command line and Python, where no
  session keeps options: it checks the options given as R does and
  returns them with the defaults for the others.

* `plotevgrid()` is in `clab_jrun()`, the command line and Python: its
  `x` maps each panel to a result of `evtst()` or `stgdid()`, or to the
  rows `evtm`, `est`, `lo` and `hi` themselves.

* `plothex()` returns its hexagons in order of row, then column, in every
  locale (it ordered them by `table()` of text keys, which sorts in the
  session's collation), and is in `clab_jrun()`, the command line and
  Python with the same order.

* `co2fp()` is in `clab_jrun()`, the command line and Python as the
  conversion it makes: the caller gives the wall and CPU seconds (and,
  where measured, the machine's CPU accounting and the RAPL counters
  before and after) and gets R's result without `value`. Timing an
  expression stays in R.

* `hwpwr()` is in `clab_jrun()`, the command line and Python. The result
  holds the rows R returns, the interval `dt` and the last reading as
  text (`snap`), which a later call takes as `since`; R keeps the last
  two as attributes.

* `plotrdg()` and `plotvln()` are in `clab_jrun()`, the command line and
  Python, returning the rows R returns (without the quartiles `plotvln()`
  keeps as an attribute). Their densities come from a kernel that
  reproduces `stats::density()`: the bandwidth rules nrd0, nrd, ucv, bcv
  and SJ, the binning, the FFT convolution and the interpolation, equal
  to R's to rounding. In R both still call `stats::density()`.

* New: `dwnexp()` builds downwind-exposure instruments from sources of
  pollution and the wind: the sector and transport-wedge rules of Bao and
  Da (2026) for smoke plumes that move from day to day, and the angle and
  distance instrument of Ayesh (2023) for fixed sources. `tailwd()` gives
  the tailwind component along a flight course, from which Hener (2022)
  builds the instrument of the landing direction for aircraft noise.

* New: `invsn()` takes `by` and `rule` to turn readings into one
  inversion per day, as Batkeyev and DeRemer (2023) require every
  three-hourly reading of a day to be inverted.

* The help pages of `windiv()`, `ivfe()`, `ppmlcf()`, `smkiv()` and
  `ccross()` cite the designs of Karamik and von Graevenitz (2022),
  Sarmiento (2020), Batkeyev and DeRemer (2023), Hener (2022), Bao and Da
  (2026) and Rahman et al. (2023), with the arguments that reproduce
  each.
* New: `dirsim()` draws compositions from a Dirichlet distribution (a
  gamma variate per part over its row sum), from R's generator, so that
  `set.seed()` and `rgamma()` give the same draws; for testing a
  compositional analysis such as `clrilr()` before real measurements exist.

* New: `evlx()` takes `true`, a true value other than the null to set the
  estimate against, as the argument of that name in the 'EValue' package.

* Fixed: the second sweep found 25 more functions (`lodimp()`, `lnroc()`,
  `ibench()`, `synexp()`, `bndxw()`, `rptvic()`, `etohlr()`, `mndla()`,
  `crchyp()`, `clrilr()`, ...) whose compiled readers asked for the length
  of `NULL`, a function or an environment; they now refuse it in words.

* Fixed, after a second sweep: a function given as `est` to the exposure
  designs, `end` of a point `datumpmap()`, a list for the data of
  `expconcidx()`, `exprspplr()` and `polbrdarea()`, `ref` and a single
  `p2` of `jntx()`, a bootstrap `seed` past an integer in the exposure
  estimators, an `x` of `sarml()` that is not numbers, an `onlsmry()`
  accumulator without its state, an infinite `cnts` of `bandx()`, a
  non-numeric `source` of `advdif2d()`, an object without a `crs` in
  `crsof()`, invalid text in `pbbnd()` and a table of no rows in
  `yoysmry()` used to fail in R's words; each is now refused in the
  package's, and `ccross()` names itself where it named `plcday()`.

* Fixed: `kscan()`'s Monte Carlo p-values used to depend on the order of
  the rows for a given seed (the times were shuffled over the rows as
  given); the shuffles now run over the incidents in the order of their
  coordinates and times. The p-values for a given seed change once; the
  clusters do not.

* Fixed: `ctabx()` of a table of zeros used to print "n = NA" and give no
  statistics; such a table has no total, and it is now refused (here and
  on the command line).

* Fixed, one function at a time: `bufas()` multiplied a factor or text
  `r` for its default `buf` before checking it (a warning or "non-numeric
  argument to binary operator"); `gtmp()` tested a function `q` for `NA`
  with a warning; the `rgnmap*()` checks turned a function into text;
  `ppmlfe()`, `ppmlcf()` and `expgmm()` converted a `maxit` past an
  integer with a warning; `yoypdf()` did not check `lag`; `dstnskch()`,
  `rsvsmp()`, `relbl()`, `grdrcd()` and `grdlvls()` took a function for
  their values or labels; `crfno2()`, `crfpm25()` and `polbrdn()` failed
  quoting a function as an outcome or pollutant. Each now refuses in its
  own words.

* Fixed: a colour index past the range of an integer (`col = 1e300`,
  `Inf`) used to pass the plots' colour check with a coercion warning, or
  fail in `col2rgb()` ("numerical color values must be >= 0"), in 12 plot
  functions; it is now refused as any bad colour ("`col` must be
  colours."). A 1 x 1 matrix for a number such as `bins` used to be
  accepted with a recycling warning; it is now refused.

* Fixed: a file path longer than a path can be (1e6 characters) used to
  reach `path.expand()` and its warning ("expanded path length ... would be
  too long") in `co2fp()`, `raplok()`, `wrtxtfb()` and the `yoy*()`
  writers; it is now refused, at most 4096 bytes.

* Fixed: text that is not valid UTF-8 given for an argument used to stop
  49 functions inside base R's string functions ("invalid multibyte
  string", "input string 1 is invalid UTF-8") or, in a title, at the
  graphics device; the shared text check now refuses it first, naming the
  argument (`bandx()`, `dstfit()`, `gplume()`, `fylbl()`, the `yoy*()`
  writers, the descriptive statistics, the plots' `main`, ...), and a
  refusal that quotes such text shows its bytes.

* Fixed: `NULL`, a function or an environment given for a vector argument
  used to reach the compiled code and fail there in R's words ("LENGTH or
  similar applied to NULL object", "cannot coerce type 'builtin' to
  vector"), in 31 functions (`calaud()`, `fairdf()`, `tmpaud()`,
  `lodimp()`, `trnmat()`, `fairsim()`, ...); it is now refused in the
  package's ("`area` must be a vector, got NULL.").
* New: `ccepnl()` estimates dynamic heterogeneous panels as Chudik and
  Pesaran (2015) set out: `ylag` lags of the outcome among the
  regressors, `lags` lags of the cross-section averages (by default the
  integer part of `T^(1/3)` when `ylag > 0`), the mean group estimator
  and its variance, and the half-panel jackknife (`bias = "jack"`) or
  recursive mean adjustment (`bias = "rma"`), with the uncorrected
  estimate alongside in `mg0`. The pooled estimator is not made for
  dynamic models (`pool`), as their theory is that of the mean group
  estimator. Calls without the new arguments give the results they gave.

* New: `ccepnl()` reports the heterogeneity-bias test of the pooled
  estimator of Juodis and Reese (2026, `hb`) and the number of factors
  the averages identify by their two regularized procedures (`nfac`),
  and offers the fixed-T variance of the pooled estimator of Brown,
  Schmidt and Wooldridge (2026, `pvar = "bsw"`). The help page sets out
  the five lessons of Juodis and Reese and what each asks of a user.

* New: `gmmpnl()` returns the residuals in levels (`resid`), for
  Pesaran's test of cross-sectional dependence and Moran's I on each
  year, and the joint test as xtabond2 prints it (`xwald`); its help page
  maps each estimate and test of Andreatta et al. (2023) to a call and has
  an example of their specification on simulated data.

* Changed: `gmmpnl(small = TRUE)` now makes the small-sample corrections
  as xtabond2 3.7.2 computes them (Roodman 2009): the conventional
  one-step variance and `sigma2` times `w / (w - k)`, the robust and
  two-step variances times `(n - 1) / (n - k) G / (G - 1)`, and t and F
  on `n - q - c` or `G - c` degrees of freedom; it used `n / (n - k)`
  throughout. The joint test counts its degrees of freedom as the
  program does, with or without `small`. `small = TRUE` is refused with
  one group, or with no more transformed rows than coefficients.

* Fixed: `xcser()` used to report r = 0 with an interval, p = 1 and
  nonzero cross-correlations for a constant series whose mean rounds
  (`rep(1, 10)`, `rep(0.1, 10)`); a correlation with a constant series is
  undefined, and `xcser()` and the command line now refuse it, naming the
  series.

* Fixed: the other statistics that divide by the spread of the data
  used to return numbers for a constant input whose mean rounds
  (`rep(0.1, n)`, `rep(1 / 3, n)`), where they refused a constant with an
  exact mean: `gearc()`, `lmorn()`, `locgi()` and `bvmor()` (Geary's C,
  local Moran, G and the bivariate Moran are undefined), `txcal()` (no
  calibration line), `rdloc()` (no Imbens-Kalyanaraman bandwidth),
  `pcdtst()` (a pair with a constant unit, now skipped) and the
  Jarque-Bera test of `fepnl()`. They now refuse it, or skip it, as for
  any constant. `locgi(star = FALSE)` gives `NaN` for the one unit whose
  value differs from all the others (its reference values are constant).
  `cmpsc()` gives no correlation to a stretch that does not vary, as
  `cor()` does, instead of one made of the rounding of its running sums.
  `ecolrg()` and `evfilt()` refuse a covariate column that does not vary
  (it cannot be told from the intercept), naming it, and `evfilt()` an
  outcome that does not vary.

* Fixed: `hsacc()` used to break ties in the predicted score at the
  cut-off of `top` by the order of the rows, so that reordering the cells
  changed the hit rate, PAI, PEI and RRI (0.571 or 0.536 at `top = 0.2` in
  one example); and `nmax` broke ties in density the same way. The new
  argument `ties` says which tied cells are hot: `"in"` (all, the
  default), `"out"` (none) or `"frac"` (each by the share of their area
  still needed); `nmax` takes the larger of equally dense cells first.
  The result no longer depends on the order of the rows; without ties at
  the cut-off it is as before.

* Fixed: the point patterns and distance functions used to run without
  end, or until the machine was stopped, on extreme radii and magnitudes:
  a radius or distance of 5e-324 (`ripk()`, `ripenv()`, `dbscn()`,
  `knxnr()` and the other functions on the shared grid of points), a
  radius of 1e9 (`ripk()`, `ripenv()`), `opteffort(W = 5e-324)`,
  `dstfit()` and `yoypx()` on counts of 1e50 and `smixm()` on coordinates
  of 1e200. The grid now starts from the smallest normal double and its
  searches stop at its edge; `opteffort()` refuses a `W` whose densities
  underflow or overflow; the goodness of fit of `dstfit()`'s discrete
  fits looks only at the jumps of the data; `yoypx()` refuses counts
  above 2^53 for its exact test; `smixm()` and `sglmx()` refuse points
  whose squared distances overflow; and the optimisers stop on an
  infinite start or search direction. Results for other inputs are
  unchanged.

* Fixed: `arimfc()` and `sarfc()` used to refuse a series that does not
  vary once differenced (a constant, a straight line) in the optimiser's
  words ("initial value in 'vmmin' is not finite"); they now say
  "`y` does not vary once differenced: the model has nothing to fit.",
  as does the command line.

* Fixed: the command line's `fylbl()` converted an `endyear` of 1e300,
  `Inf` or `NaN` to a whole number with undefined behaviour; like R it now
  gives `NA` for a year outside the range of an integer. `ousim()`
  counted `days / dt` steps as an integer before checking them (1e16 at
  `dt = 1e-9`), and the periodic baseline of the Hawkes functions its
  integration steps; both are now counted in double first.

* Fixed: the command line and `clab_jrun()` used to copy text that is not
  valid UTF-8 from an argument or a CSV cell into their output, which was
  then not valid JSON; they now refuse it, naming the argument, in the
  words of R ("`group` holds text that is not valid UTF-8."). A NUL byte
  in a CSV cell used to end the cell's number silently (`4<NUL>0` was
  read as 4); it is now refused with its row and column.

* Fixed: a function that calls another (`runfx()` calls `exprun()` and an
  estimator, `plmfld()` `gplume()`, `polbrdn()` `crfno2()`, `rsvsmp()`
  `rsvidx()`, `xfrvrf()` `vrfmrg()`, `yoywr()` the writers, and the
  exposure functions `ppmlfe()`) used to pass on the inner function's
  refusals and warnings under the inner name ("gplume: ..."); they now
  carry the name of the function called (`ccross()`, `dosebin()`,
  `hetfx()`, `jntx()`, `lftab()`, `plcday()`, `plmfld()`, `polbrdn()`,
  `prechk()`, `rsvsmp()`, `runfx()`, `typefx()`, `xfrvrf()`, `yoyhtml()`,
  `yoywr()`).
* Fixed: the choice arguments of the exposure functions (`ppmlfe()`,
  `ppmlcf()`, `ivfe()`, `ccross()`, `expgmm()`, `hetfx()`, `invsn()`,
  `jntx()`, `plcday()`, `rdloc()`, `smkiv()`, `typefx()`, `updwn()`,
  `windiv()`) used to be refused by `match.arg()` in base R's words
  ("'arg' should be one of ..."); they are now refused like every other
  choice ("ppmlfe: `se` must be one of ...").

* Fixed: `hwpwr()` used to report as measured a RAPL counter that fell
  by more than a wrap could explain (a reset read as 257 kJ) or read above
  its range, and read counter text that is not a whole decimal number
  ("0x10", "1.5e3", an embedded NUL, 26 digits); such rows are now not
  measured, with the reason in `why`. A FIFO in place of a counter used to
  block the call for ever; files are now read only when regular, opened
  without blocking. A library that is not NVML, or an `nvidia-smi` that
  fails, no longer adds a GPU row. `dt` is at most 604800 (a week) and
  `pid` at most 2147483647 (`pid = 2^31` used to be accepted with a
  warning).

* Environmental crime, with offences against the environment as the
  outcome: the LMDI decomposition of a count over the factors of a Kaya
  identity, over time and against a reference region (`lmdix()`, after
  Duan et al. 2026), the two-way fixed-effects regression of a count on
  drivers with an interaction, standardised coefficients, variance
  inflation factors, simple slopes and the Jarque-Bera test of the
  residuals (`fepnl()`), and Pesaran's CD test of cross-sectional
  dependence (`pcdtst()`, as Andreatta et al. 2023 use it). They are in
  `clab_jrun()`, the command line and Python.

* Dynamic and heterogeneous panels, for the models Andreatta et al. (2023)
  fit to waste crime: difference and system GMM (`gmmpnl()`, Arellano and
  Bond 1991, Arellano and Bover 1995, Blundell and Bond 1998, with the
  instruments of Roodman 2009), one-step and two-step, first differences
  or forward orthogonal deviations, GMM-style instruments with lag limits
  and collapsing, IV-style instruments, the robust or conventional
  one-step variance and Windmeijer's corrected or the uncorrected two-step
  variance, small-sample t and F statistics, a joint Wald or F test, the
  Arellano-Bond tests, the Sargan, Hansen and difference-in-Hansen tests
  and the instrument count; Wooldridge's test of serial correlation
  (`wsctst()`, Drukker 2003); and the common correlated effects mean group
  and pooled estimators (`ccepnl()`, Pesaran 2006), on balanced and
  unbalanced panels, with a choice of the averages taken out. They are in
  `clab_jrun()`, the command line and Python.
* Panels with a repeated unit and period are found by sorting the rows,
  not with a grid of units by periods, which for many units and periods
  could ask for more memory than the rows hold. `pcdtst()` and `lmdix()`
  refuse a panel whose units by periods pass 2^26 cells.

* Spatial epidemiology and surveillance: Bayes outbreak detection
  (`bayout()`), exposure within a buffer and by kernel (`bufexp()`,
  `kerexp()`), the Leroux and BYM2 area structures (`lerprc()`,
  `bym2sc()`), ecological Poisson, negative binomial and zero-inflated
  regression with expected counts (`ecolrg()`), eigenvector spatial
  filtering (`evfilt()`), the metapopulation SIR model (`spsir()`), Weibull
  frailty and promotion-time cure models with CAR area effects (`sfrail()`,
  `sfcure()`) and the spatial encoding of graph attention (`spbias()`).
  They are in `clab_jrun()`, the command line and Python.

* Geostatistics on planar coordinates: the empirical semivariogram
  (`vgemp()`), variogram models fitted by maximum likelihood (`vgfit()`),
  ordinary kriging (`okrge()`), block kriging (`bkrig()`) and its
  leave-one-out cross-validation (`krigcv()`), polynomial trend surfaces
  (`trsurf()`), the spatial linear model by profile likelihood (`sglmx()`)
  and the spatial linear mixed model by REML (`smixm()`), with a spatial
  data object (`spobj()`, `spxy()`, `crsof()`), the specification of a web
  map (`lfspec()`) and density-weighted sampling of points (`spsamp()`).
  They are in `clab_jrun()`, the command line and Python.

* The energy and carbon footprint of a computation (`co2fp()`), modelled
  from running time, utilisation and memory (Green Algorithms) or measured
  from the CPU's RAPL counters (`raplok()`), with the measured energy of
  other devices added from a table; the carbon intensity of electricity by
  location, offline (`co2int()`, `co2inttab()`), the location found
  without network access (`dtctloc()`) and everyday equivalents
  (`fpeqv()`). All but `co2fp()`, which evaluates an R expression, are in
  `clab_jrun()`, the command line and Python.

* Air pollution, for environmental criminology: the exposure-response
  estimate with a bootstrap sensitivity interval (`exprspplr()`), by the
  cross-fitted partially linear regression in the compiled core or any
  estimator with its interface, and the whole pollution-to-health pipeline
  from demo data, a (compressed) CSV or scalars to a report and its text
  (`vrfpol()`, `polrpttxt()`). `exprspplr()` with its built-in estimator is
  in `clab_jrun()`, the command line and Python; the other two read files
  and write text in R.

* Scaling and ideal points, for judges, courts, legislators and survey
  respondents: vote agreement (`votagr()`), Aldrich-McKelvey scaling
  (`amscal()`) and its Bayesian form (`bamscl()`), Basic Space scaling
  (`bbscal()`), classical, SMACOF and nonmetric MDS (`cmdsf()`, `smacf()`,
  `nmdsf()`, `dblctr()`, `mdsfit()`), Bayesian MDS and unfolding
  (`bmdsf()`, `bunfld()`), unfolding by least squares and by SMACOF
  (`mlsmu()`, `smunf()`, `unfstr()`, `ipsrec()`), INDSCAL (`indscl()`),
  Procrustes rotation (`procrs()`), NOMINATE utilities, likelihood, fit
  and bootstrap (`nomutl()`, `nompr()`, `nomll()`, `dwnom()`, `nombt()`),
  alpha-NOMINATE (`anomin()`), Optimal Classification and its ordered
  form (`optcls()`, `ordocl()`), binary IRT by EM (`irtem()`) and by Gibbs
  sampling (`cjrirt()`, `irtlik()`, `irtpst()`), ordinal and dynamic IRT
  (`ordirt()`, `dynirt()`), anchoring vignettes (`anchv()`), the bootstrap
  of scaling positions (`scboot()`), normal vectors and cutting lines
  (`nrmvec()`, `cutlns()`) and Wordfish (`wrdfsh()`). The samplers draw
  from R's generator seeded as `set.seed()` (`cjrirt()` from Philox), so a
  chain is the same in every front end. All are in `clab_jrun()`, the
  command line and Python.

* Counts over time and across categories: ARIMA and seasonal ARIMA
  forecasts by exact maximum likelihood (`arimfc()`, `sarfc()`), the
  Ornstein-Uhlenbeck fit, simulation and Fokker-Planck density of a count
  series (`ousim()`, `fpdens()`), Lotka-Volterra rates of a yearly series
  (`lvfit()`) and the correlation of counts per unit between categories
  (`catcor()`). They are in `clab_jrun()`, the command line and Python.

* Regression discontinuity: local-linear sharp and fuzzy estimates with
  the bandwidth of Imbens and Kalyanaraman (2012) or one given, three
  kernels, conventional and robust intervals (Calonico, Cattaneo and
  Titiunik 2014, with the pilot bandwidth equal to the main one),
  heteroskedasticity-robust or clustered errors and placebo cutoffs
  (`rdloc()`); the density test of McCrary (2008) (`rddens()`).

* Exposure and crime, checked against the papers: the Newey-West lag
  of `se = "nw"` when `lag` is not given is now the automatic bandwidth
  of Newey and West (1994, equation 2.2) instead of their lag-selection
  parameter `n`, so results with that default change; two-way clustering
  gains the factor of each component (`twc = "each"`) and the
  eigenvalue fix (`psd = "eig"`) of Cameron, Gelbach and Miller (2011)
  in `ppmlfe()` and `ivfe()`. The first-stage table of `ivfe()` gains the
  Kleibergen-Paap rk statistic for several endogenous columns (it was
  `NA`) and the effective degrees of freedom and 5% critical value of
  the simplified test of Montiel Olea and Pflueger (2013); the effective
  F itself was checked against their equation 3 and is unchanged.

* `ppmlfe()` and `ivfe()` absorb any number of sets of fixed effects
  (a third set such as hour or region-by-year, which several of the
  papers use); one or two sets are computed exactly as before.

* `ccross()` gains single-day lags (`lagv`, `lag`) and natural-spline
  controls (`ns`), the specification of the nationwide case-crossover
  study of Han et al. (2026), with a worked example of its main table.

* `smkiv(cf = "childs")`: smoke PM2.5 as Childs et al. (2022) define it,
  PM2.5 above the three-year median of the same place and month on
  smoke-free days, on smoke days.

* `ppmlfe(se = "quasi")`, the quasi-Poisson variance with the Pearson
  dispersion, the overdispersion adjustment Armstrong, Gasparrini and
  Tobias (2014) give for the case-crossover design of `ccross()`.

* `ivfe()` reports the Sargan statistic and the cluster-robust Hansen J
  of the two-step GMM estimate when there are more instruments than
  endogenous columns (`overid`, `gmm`), the cross-definition
  over-identification test the review of Yamada (2026) recommends.

* `ppmlfe(sep = "sst")`: the check of Santos Silva and Tenreyro (2010)
  for Poisson estimates that do not exist, dropping the regressors that
  are collinear on the rows with positive counts. The iterative
  rectifier of Correia, Guimaraes and Zylkin, which drops rows instead,
  is not built: their 2020 paper names it without stating it.

* Exposure and crime, a count with an endogenous exposure by the
  exponential-mean generalised method of moments of Mullahy (1997) and
  Windmeijer and Santos Silva (1997): multiplicative or additive
  residual, one-step, two-step or iterated weighting, cluster-robust
  errors and Hansen's J (`expgmm()`).

* Exposure and crime, the designs of 2019-2026: smoke instruments from
  plume data, raw or adjusted to the plumes that reached the ground
  (`smkiv()`); runs of exposed days and the effect by day in the run
  (`exprun()`, `runfx()`); effects by a moderator with the test of
  equality (`hetfx()`); the time-stratified case-crossover design
  (`ccross()`); two exposures together as a surface (`jntx()`); exposure
  in bins with tests of non-linearity (`dosebin()`); and effects by
  offence type with a multiple-testing correction (`typefx()`). Every
  choice the papers leave open is an option of the estimators, with the
  standard errors of each kind, the three first-stage statistics and the
  plug-in of the count instrument.

* Exposure and crime, for counts by place and day against a daily
  exposure: Poisson pseudo-maximum likelihood with one or two sets of
  absorbed fixed effects, an offset and one- or two-way clustered errors
  (`ppmlfe()`); two-stage least squares with the same effects and the
  first-stage Wald statistic (`ivfe()`); the Poisson control function with
  a cluster bootstrap (`ppmlcf()`); wind-direction and inversion
  instruments (`windiv()`, `invsn()`); the irrelevant-day placebo
  (`plcday()`); the upwind and downwind design around a road or a point
  source (`updwn()`); and the join of exposure to counts (`expjn()`). The
  methods are those of Bondy, Roth and Sager (2020) and of Herrnstadt,
  Heyes, Muehlegger and Saberian (2019, 2021).

* Air pollution, for environmental criminology: Pasquill-Gifford
  dispersion coefficients (`pgsig()`), Briggs plume rise (`briggs()`), the
  Gaussian plume and puff (`gplume()`, `gspuff()`, the field
  `plmfld()`), an advection-diffusion grid (`advdif2d()`) and a particle
  walk (`lagrprt()`); the concentration-response functions of PM2.5 and
  NO2 (`crfpm25()`, `crfno2()`), the attributable fraction
  (`attrfrac()`), deaths avoided (`mortdisp()`), the burden of a
  population and of areas (`polbrdn()`, `polbrdarea()`) and the
  concentration index of exposure by income (`expconcidx()`).

* Describing a table before trusting it: inline histograms
  (`inlhist()`), frequency tables (`freqtbl()`), pairwise correlations
  (`cortbl()`, `topcors()`), empty and constant rows and columns
  (`drpempty()`, `drpconst()`), clean column names (`clncolnms()`),
  multivariate outliers (`mahaoutl()`), the shape of missingness
  (`missruns()`, `missgrid()`, `misspat()`, `misssmry()`) and Little's test
  of missing completely at random (`mcartst()`), duplicated rows
  (`duprows()`), a column report (`profcols()`), candidate distributions
  fitted by maximum likelihood and ranked by AIC (`dstfit()`) and
  Cochran's Q for distinct groups before pooling (`dstntst()`). They are
  in `clab_jrun()`, the command line and Python.

* One-pass summaries of data larger than memory: moment accumulators
  (`onlsmry()`, `smryupd()`, `smrymrg()`, `smrystat()`), reservoir sampling
  (`rsvidx()`, `rsvsmp()`) and the HyperLogLog distinct count
  (`dstnskch()`, `dstncnt()`, `skchmrg()`).

* Search theory for search and rescue planning: lateral range curves
  (`ltrlrng()`) and sweep width (`swpwdth()`), the probability of detection
  from coverage (`srchpod()`) and of a parallel sweep with navigation error
  (`pswppod()`), the probability of success with Bayesian updating of
  containment (`srchscs()`), datum probability maps (`datumpmap()`), the
  optimal allocation of effort (`opteffort()`, `optcrcsrch()`), Frost's
  optimal search squares (`optsqr()`) and the list of his book's figures
  (`srchfig()`). They are in `clab_jrun()`, the command line and Python.
* First Nations, Metis and Inuit health and wellbeing, a core part of the
  package beside the justice system it concerns: Statistics Canada's
  published tables read with their quality flags (`sctidy()`, `scests()`,
  `qltyrcd()`, `limfrom()`), compared group by group and as whole rated
  scales (`indgap()`, `ordcmp()`), standardised by age (`stdprev()`),
  wellness beside distress (`welprof()`, `dualcont()`, `holbal()`,
  `frmaud()`), the First Nations Mental Wellness Continuum (`mhcont()`),
  healing as change (`healchg()`), connection across places
  (`connassoc()`), community-weighted indices (`commidx()`) and the
  Community Well-Being index (`cwbscr()`); the generations of residential
  school attendance and their dose-response (`genexpo()`, `exptrnd()`);
  the 2021 Census Indigenous Population Profile (`cnsprof()`, `profshr()`,
  `cnsstd()`, `cnscwb()`); and the Calls to Action (`trcind()`,
  `trchlthgap()`, `trccmpl()`, `trcprog()`). They are in `clab_jrun()`, the
  command line and Python.

* `polex()` makes the simulated policy panel the examples of the policy
  functions use: twenty districts, 36 months and a staged reform with a
  known effect. The same seed gives the same panel in every front end.
* `thrtst()`: when the quasi-Newton steps of the maximum-likelihood fit
  end short of the gradient criterion, up to 20 Newton steps now follow,
  so the fit converges on flat likelihoods on every architecture (arm64
  builds reported the self-test's fit as not converged). Fits that met
  the criterion give the same results as before.

* `thrtst()` gains `fit = "bayes"`: the threshold test fitted by Bayesian
  inference, as Simoiu, Corbett-Davies and Goel (2017) and Pierson,
  Corbett-Davies and Goel (2017) fit it, with a threshold for every group
  and area and the hierarchical priors that Simoiu et al. state, on the
  discriminant distributions. The posterior is sampled in C++ by the
  No-U-Turn sampler (`chains`, `iter`, `seed`), starting near the
  maximum-likelihood fit. `thr` gives the posterior means and intervals
  of the thresholds and of their differences from the reference group,
  and `post` the summaries, effective sample size and R-hat of every
  parameter. Chains with an R-hat of 1.05 or more, or fewer than 100
  effective draws of a threshold, are refused. The fit agrees with Stan's
  sampler on the same model. `fit = "ml"` stays the default. The two fits
  answer under different models and should be reported as such, not as
  checks of one another: with a threshold for every group and area under
  a hierarchical prior, on data made with one threshold per group the
  Bayesian fit can report thresholds far from the maximum-likelihood ones
  and from the generating values (three groups in ten areas, 20,000 stops
  a cell: about 0.13, 0.07 and 0.14 against 0.27, 0.20 and 0.31), which
  Stan reproduces on the same model.

* `thrtst()` gains `rank = "joint"`, the rank-order test of Anwar and
  Fang (2006) as the paper states it: the ranking of all officer groups
  is tested at once by their one-sided z statistics, and the null of no
  relative prejudice is rejected when the significant orders form a
  cycle; the common ranking is given otherwise. The pairwise form stays
  the default. With `off`, the new element `mono` holds their test of
  monolithic behaviour (Pearson's chi-square of equal rates across
  officer groups, by group of the people stopped and measure).

* `ibench()` gains the two methods of Ridgeway and MacDonald (2009) it
  left out: `dr = TRUE` gives the doubly robust benchmark (the
  propensity-weighted logistic regression of the outcome on the officer's
  indicator and the circumstances, with a sandwich standard error), and
  `local = TRUE` Efron's local false discovery rate of each officer, with
  the empirical null fitted as the package locfdr does by default (column
  `lfdr`, element `null`), flagging by it. Benjamini and Hochberg's `q`
  stays. The defaults give the results as before. The tests check the
  doubly robust benchmark against `glm()` and a sandwich by hand, and the
  local rate against locfdr (reference values from `dev/ref/lfdr.R`).

* `lnsize()` follows Malpass (1981) for `esize`: members nobody chose
  (null foils) leave the nominal size and the chance expectation, as
  Tredoux (1996, Table 4.2) shows. Before, `esize` did not adjust; that
  value is now `esizeb`. The two differ only when a member was never
  chosen. The interval of Tredoux's `E'` (`etrlo`, `etrhi`) now takes the
  ends of the interval of `I = 1 - sum p^2` to `1 / (1 - I)`, as Tredoux
  does, in place of a symmetric delta-method interval around `E'`; the
  ends move for every lineup (for 40 mock witnesses choosing 15, 5, 6, 4,
  7 and 3 of six, E' 4.44: from 2.98 to 5.90 before, 3.35 to 6.62 now). A choice
  of 0 records a mock witness who chose no one: it counts for `psusp` and
  `fsize` (as in Wells, Leippe and Ostrom 1979) and not for the spread of
  choices; `nm` is the number who chose a member. The tests reproduce the
  worked examples of Tredoux (1996, 1999).

* `lnroc()` continues a curve that ends before `far` by the straight line
  to the point where every lineup counts as an identification, as pROC
  does, in place of a flat line; this changes `pauc` only when `far` is
  given above a procedure's largest false alarm rate (the default never
  is). The tests reproduce Table 3 and the partial areas of Mickes, Flowe
  and Wixted (2012).

* `wdqx()` now names the crimes prevented in the treated area and the
  buffer `tne`, the total net effect of Guerette (2009, Appendix B,
  equation 4), in place of `net` (same value), and adds the gross effect
  `ge` and the net effect `ne` of the same appendix. Its help page gives
  the formulas and Guerette's worked example; the tests reproduce it.

* Every result that holds something to draw has a `plot()` method on base
  graphics: estimates with intervals as dot-and-interval charts, effects by
  time since adoption, paths, survival and ROC curves, Lorenz curves, bars,
  matrices with their numbers in the cells, and surfaces. Each returns,
  invisibly, the data it drew. `blpal()` gives the colour-blind-safe
  palettes they share, `plotopt()` their session look and `yoypals()`
  the palettes of change tables; series differ in
  line type or point shape as well as colour. The general plots come with
  them: Lorenz and funnel plots (`plotlrnz()`, `plotfnl()`), distributions
  by group (`plotrdg()`, `plotcdfs()`, `plotvln()`, `plothex()`,
  `plotgain()`), series over time (`plotcal()`, `plothrz()`, `plotmult()`,
  `plotbump()`) and event studies in a grid (`plotevgrid()`).
  `plot()` of an `hwkres` result now returns the quantiles it drew.
* Change tables from `yoypx()` (or any table with its columns) can be
  labelled (`yoylbl()`, `fylbl()`), summarised (`yoysmry()`) and written:
  CSV, TSV and JSON with every digit and every code as it is (`yoycsv()`,
  `yoytsv()`, `yoyjs()`), Markdown (`yoymd()`), a self-contained HTML page
  (`yoyhtml()`) and a PDF written by the package itself (`yoypdf()`), the
  format picked from the file name by `yoywr()`.

* The drift battery: the two-sample Kolmogorov-Smirnov test (`drftks()`),
  the population stability index and Jensen-Shannon divergence
  (`drftpsi()`) and the chi-square tests of a categorical column against a
  known distribution and of two samples' homogeneity (`drftchisq()`,
  `drfthomog()`), each result keeping its inputs in the attribute
  `lakner`. They check a new release of a published table against the last
  one before year-over-year change is computed from it. They are in
  `clab_jrun()`, the command line and Python (there without the Monte
  Carlo branch of `drfthomog()`).

* `rgnmapintg()`, `rgnmapcmp()` and `rgnmapalt()` check a crosswalk from
  units to regions: each unit once and with a known region, a rebuilt map
  against the published one, and the map against an assignment by another
  route with the known disagreements named. They are in `clab_jrun()`, the
  command line and Python.

* The categorical-integrity checks: codes decoded by an explicit
  dictionary or by their value labels (`deccodes()`, `declbld()`), recodes
  and relabels by name with every value mapped (`grdrcd()`, `relbl()`),
  declared levels (`grdlvls()`), a 0/1 indicator that refuses a factor
  (`grdbin()`), a recode proved by its cross-tabulation (`vrfrcd()`),
  counts checked against those a release published, with the label
  permutation that would explain a mismatch (`vrfmrg()`), and a variable
  checked after it crossed from one program to another (`xfrvrf()`). They
  are in `clab_jrun()`, the command line and Python.

* `shrs()` gives each count's share of a total in percent with Wilson's
  score interval, from a vector or from a table summed by groups. It is in
  `clab_jrun()`, the command line and Python.

* `toascii()`, `asciifb()` and `wrtxtfb()` transliterate text to plain
  ASCII by a fixed table, keep valid UTF-8 as it is, and write text as
  UTF-8 bytes whatever the locale. They are in `clab_jrun()`, the command
  line and Python.

* The summary statistics of the compiled core are exported: the mean,
  variance, correlation and standard deviation (`cmn()`, `cvari()`,
  `ccorr()`, `csdev()`, `cdistm()`), the normal density (`cnormpdf()`,
  `cnormlpdf()`), moments (`cmom()`), quantiles and robust spread
  (`cqnt()`, `cmed()`, `cmad()`, `ciqr()`, `ctkyfnc()`), trimmed and
  winsorized means (`ctrmmean()`, `cwinsmean()`), weighted moments
  (`cwtd()`), ranks (`cspr()`, `cmidrk()`), the covariance matrix
  (`ccovm()`), bootstrap means (`cbootmean()`), trimmed inverse-
  probability weights (`cipw()`), the incomplete gamma function
  (`cgamcdf()`) and splitmix64 uniforms (`cunif()`). They are in
  `clab_jrun()`, the command line and Python.

* `prechk()` checks a table against a preset before `clprep()` reads it,
  and names in one error every column the preset reads (directly, in a
  filter or in a recipe) that the table lacks, every label outside the
  codes the publisher lists for its column, and every label a mapping
  recode without a default does not cover. The codes of the Ontario
  segregation, restrictive confinement and deaths-in-custody tables are
  `inst/presets/codes.tsv`, from the ministry's data dictionary; a preset
  of one's own takes its codes as a data frame. The presets describe the
  layouts of published tables, and crimlab ships none of the tables.

* The OTIS treatment-outcome presets of `cbatt()` (`otis-a01-cbatt`,
  `otis-a01n-cbatt`, `otis-pa-cbatt`, `otis-pb-cbatt`, `otis-pc-cbatt`)
  are tested against rmorie's pair builders on a simulated placement
  table, with the values stored in `tests/testthat/ref/preset.txt` by
  `dev/ref/preset.R`.

* `hwpwr()` reads the power and energy of the machine's devices over an
  interval from the operating system's own counters: RAPL on Linux, the
  NVIDIA management library, the amdgpu and Intel GPU hwmon files, the
  IOReport energy counters of Apple silicon, the Energy Meter Interface
  on Windows, and the battery while it discharges. Each figure says how
  it was obtained; an estimate carries its range, and a device that
  cannot be measured is listed with the reason.

* `mtrbnd()` bounds the mean outcome under each of ordered sentences and
  the effect of each against the mildest, under the monotone treatment
  response and selection assumptions of Manski and Pepper. It is in
  `clab_jrun()`, the command line and Python.

* `hsacc()` scores a hot-spot map on later crimes: the hit rate, the
  predictive accuracy index, the predictive efficiency index and the
  recapture rate index. It is in `clab_jrun()`, the command line and
  Python.

* `apcie()` fits the age-period-cohort model of a table of counts by age
  and period with the intrinsic estimator of Yang, Fu and Land: the
  Poisson log-linear fit projected off the null vector of the design.
  It is in `clab_jrun()`, the command line and Python.

* `spcure()` fits the split-population (cure) model of Schmidt and Witte:
  a logistic share who would ever fail and a Weibull, exponential or
  lognormal time to failure among them, with covariates in both parts.
  It is in `clab_jrun()`, the command line and Python.

* `lnkpr()` links crimes into series: the distance, time and Jaccard
  similarity of behaviour of every pair, a logistic regression of being
  linked on them over the solved pairs, the area under the ROC curve of
  each, and the chance of a link for every pair. It is in `clab_jrun()`,
  the command line and Python.

* `prohs()` forecasts the risk of each grid cell from recent nearby
  incidents, weighted down with distance and age (Bowers, Johnson and
  Pease), and scores its hot cells on later incidents. It is in
  `clab_jrun()`, the command line and Python.

* `wdqx()` gives the weighted displacement quotient of Bowers and Johnson,
  with the counts expected after and the crimes prevented in the treated
  area and its buffer. It is in `clab_jrun()`, the command line and
  Python.

* `rptvic()` measures repeat victimisation: incidence, prevalence and
  concentration, the victims by number of victimisations against a
  Poisson count over the targets at risk, and the time between a
  target's successive incidents. It is in `clab_jrun()`, the command
  line and Python.

* `lnroc()` gives the ROC curve of eyewitness identifications by
  confidence, the partial area with a bootstrap comparison of procedures,
  the diagnosticity ratio and the confidence-accuracy characteristic;
  `lnsize()` the fairness of a lineup from mock witnesses (functional
  size, Malpass's and Tredoux's effective sizes, the suspect's share).
  Both are in `clab_jrun()`, the command line and Python.

* `ibench()` benchmarks each officer's stops against other officers'
  stops in the same circumstances (Ridgeway and MacDonald), by propensity
  weights, with the Benjamini-Hochberg false discovery rate over the
  officers. It is in `clab_jrun()`, the command line and Python.

* `thrtst()` tests discrimination in searches three ways: the hit-rate
  test of Knowles, Persico and Todd, the threshold test with the
  discriminant distributions of Pierson, Corbett-Davies and Goel fitted
  by maximum likelihood, and a pairwise form of the rank-order test of
  Anwar and Fang. It is in `clab_jrun()`, the command line and Python.

* `gmixlr()` now follows EuroForMix's model where they differed: only the
  alleles seen are kept apart, every other allele (those of the reference
  profiles not seen among them) is one allele that does not stutter, sized
  as the largest size listed at the marker, and the degradation slope may
  exceed 1. Its log likelihood equals EuroForMix's to 1e-6 at given
  parameters on eight cases and at the maximum on seven (one and two
  contributors, drop-out, drop-in, stutter, degradation, coancestry;
  `dev/ref/efm.R`);
  the rest is EuroForMix's single precision.
* `cmpsc()` now reproduces cmpsR exactly: a segment widened at the ends
  of the signature is cut, not moved; a lag counts when enough points
  overlap (so a segment may hang over the end of the reference); every
  confirmed peak is a candidate; ties are broken as cmpsR breaks them; and
  lags may fall on half points. `full = TRUE` keeps the earlier rule of
  lags wholly inside the reference. The scores, segments, lags and refusals
  equal cmpsR's on 16 stored cases, edge cases among them.
* `pedlr()`: the likelihood ratio of two pedigrees, given as tables of
  people and their parents, from DNA profiles at autosomal markers. The
  people not typed are summed over exactly by variable elimination on the
  alleles each received from each parent, so loops, inbreeding and any
  number of generations are allowed; mutation follows the equal or the
  stepwise model, with one rate or separate maternal and paternal rates.
  It equals `kinidx()` and `patidx()` where they apply, the enumeration of
  every allele in its tests, and pedprobr's likelihoods.
* `gmixlr()`: the likelihood ratio of two accounts of a DNA mixture by the
  gamma model of peak heights, with drop-out, drop-in, backward stutter,
  degradation and coancestry, the unknown contributors summed over and the
  parameters fitted by maximum likelihood under each account, or fixed
  with `par`. Its likelihood is checked against enumeration in R; the help
  page states the model's assumptions and limits.
* `cmpsc()`: the congruent matching profile segments score of two bullet
  signatures, with the multi-segment-length peaks of Chen et al. (2019).
  It follows the algorithm written out in R in its tests and is within one
  segment of cmpsR on simulated signatures.
* `clprep()` reaches the command line, the Python package and the
  bindings through `clab_jrun()`, which does its recodes, filters and
  recipes in C++ with R's messages and returns the generic's arguments as
  `args` with the attributes of the R result beside them. Its preset table
  is `cli/pre.inc`, made from `inst/presets/presets.tsv` by `cli/pre.sh`;
  CI fails when the two differ. `cli/test.tsv` holds a case for every
  preset. Outside R, text is ordered by its bytes (R's order under
  `LC_COLLATE=C`) and lowered as ASCII.
* `clab_jrun()`, and with it the command-line program, the Python package
  and the bindings, now reaches every exported function:
  counts, rates and tables, the audits of risk scores, judges' leniency,
  long stays, published figures and bands, stocks and flows, the Hawkes
  fits and other events in time, difference in differences, synthetic
  controls and the event study, surveys, few-cluster inference, the
  double machine learner and the super learner, and the spatial methods
  of near repeats, scans, densities, Moran's I and spatial regression. A
  warning R gives is returned in the result as `warnings`, an array of its
  messages. `lrpost()`, `impfpr()` and `lrcal()` there now give one value
  as a number, not an array of one, as R does. `src/dst.cpp` adds the
  chi-square (central and noncentral), t, F, gamma, beta, Poisson and
  binomial distribution functions, upper tails and quantiles the front
  ends need; the R functions keep calling R's own. The test of Ctrl-C in
  a running session measures the stop from the signal and sets the
  signal only once the call's handler is in place, so a loaded machine no
  longer fails it.
* The command-line program, the Python package and the Go, Rust, Java and
  JavaScript bindings now reach most of the package's functions through
  `clab_jrun()`: fairness audits, trace evidence, aoristic weights and
  digit tests, sensitivity bounds, trajectory models, networks, the
  reporting checks, the spatial and hot-spot methods, spatial weights,
  toxicology, deterrence and severity, concentration and the dark figure,
  crosswalks and taphonomy, the learners, risk terrain and the life
  course, covariate balance, judge designs and Hawkes processes. Each
  reads its arguments with the R function's checks and words and returns
  R's result as JSON; `cli/test.tsv` holds the values of the R package for
  them (`dev/clref.R`). The normal quantile and distribution function the
  front ends need are kernels of their own (`src/dst.cpp`).
* The tests no longer call other packages: the values of circular,
  DescTools, geosphere, geepack, igraph, MASS, sandwich, sf, spatstat,
  spdep, rmorie and rmoriebricklayer they compare with are stored in
  `tests/testthat/ref`, with the versions that computed them, by the
  scripts in `dev/ref`. Suggests now holds only testthat, knitr and
  rmarkdown. The `ref` workflow runs those scripts again every week on
  Linux (x86-64 and arm64), macOS (arm64 and Intel) and Windows, and
  fails when a value its platform's packages give differs from the
  stored one by more than the tolerance of the test that uses it.
* `suplrn()` stacks the learners of the package (least squares or
  logistic regression, the lasso, ridge and elastic net, the random
  forest) into a super learner: K-fold cross-fitted predictions, weights by
  non-negative least squares that sum to one, or the discrete super learner
  (the learner of least cross-validated risk), for numeric and 0/1
  outcomes, the folds drawn from a seed and returned. With one learner the
  stack is that learner. `dmlx()` takes it as a nuisance learner
  (`lrn = "stack"`, with `stk`, `nstk` and `disc`).
* Ctrl-C (Esc in the R GUIs) now stops every long computation: the
  resampling, permutation, bootstrap and simulation loops, the iterative
  fits (elastic net, logistic and count regressions, EM, quasi-Newton,
  Weber point, k-means, power iterations, Louvain, synthetic control
  weights, fixed-effects sweeps) and the scans over pairs check for an
  interrupt every few milliseconds of work, free what they hold and end
  with R's interrupt, never with a half-made result. In 'Python' the same
  loops raise `KeyboardInterrupt`.
* `etasf()` fits the temporal ETAS model (marks and the modified Omori
  law), `mvhwk()` the multivariate Hawkes process of several types that
  excite each other, with its excitation matrix and spectral radius, and
  `sthwk()` the space-time Hawkes process of crime, with intensities at
  query points; `hwksim()` simulates from a space-time fit.
  `proof("etasf")` names the theorem that ties the Omori kernel to the
  Lomax kernel of `hwkfit()`.
* `sdbsim()` simulates the Short-D'Orsogna-Brantingham model of
  burglary hot spots, as agents on a lattice or as its reaction-diffusion
  limit with suppressed sites, and `sdbthr()` gives the linear stability
  of its uniform state: the threshold beyond which hot spots form and the
  fastest-growing wavelength. `proof("sdbthr")` and `proof("sdbsim")` name
  the Lean theorems on the steady state, the trace and the decay.
* `hypq()` solves Larson's hypercube queue of patrol units exactly
  (workloads, dispatch shares, lost or waiting calls, travel times) and
  `sqrtlw()` gives Kolesar and Blum's square-root law of response distance;
  `proof("sqrtlw")` names its scaling theorems.
* `kmrec()` gives Kaplan-Meier curves by group with Greenwood's standard
  errors, `survfit()`'s log-scale intervals and medians, and the log-rank
  test; with `unit` it takes event times and studies the gaps between a
  unit's events (repeat incidents at a place), the last one censored.
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
* `crsev()` builds a severity-weighted crime index from counts by offence
  (the Crime Severity Index of Statistics Canada: weighted counts per
  population against a base, the weights given or made from incarceration
  rates and mean sentences); `zcomp()` ranks units on a composite of
  standardised counts across categories.
* `cconc()` measures how concentrated a total is across units: the Gini
  coefficient and its exact split into the share of units at zero and the
  Gini of the rest, the Lorenz curve, top shares, dispersion against the
  Poisson, the bound chance alone allows and a simulated Poisson null, the
  Gini within groups, the urban-scaling exponent, and for the places of
  events in order the growth of distinct places. `hilltl()` estimates a
  power-law tail index: the discrete maximum likelihood estimate through the
  Hurwitz zeta function, or the Hill estimators, with the Kolmogorov-Smirnov
  distance of the fit. `proof("cconc")` names the Lean theorems on the split
  at zero.
* `crcap()` estimates a hidden population from two lists (Lincoln-Petersen
  with its bounds over a box of list dependence, Chapman's estimate with
  Seber's standard error) or three (the never-listed cell under no three-way
  interaction, its sensitivity to that interaction, and the pairwise
  estimates); `dkfig()` bounds the dark figure of crime: the true
  victimisation rate under limits on a survey's false positives and
  negatives with the breakdown point of a conclusion, the total when a
  survey covers only some offences, and the offences behind incidents
  counted under the hierarchy rule. `proof("crcap")` and `proof("dkfig")`
  name the Lean theorems.
* `bandx()` reads band labels ("1-5", "Over 60", "Under 10") into bounds,
  places a value in each band, expands a banded frequency table and shows
  how much a statistic depends on the cap assumed for an open band;
  `pbbnd()` gives the bounds a rounded or suppressed published cell allows
  and the envelope of the change between two published columns; `yoypx()`
  gives year-over-year (or any period) change within groups with the exact
  interval of the percent change, exact p-values adjusted for multiplicity,
  and the envelope that rounding allows. `proof("pbbnd")` names the Lean
  theorems on the envelope.
* `jdgiv()` gives the judge-leniency instrument's first stage, intention-to-
  treat effect and Wald ratio with the compliers' effect over a grid of
  defier shares and effects, or the exact decomposition by type on a known
  population; `jmono()` gives each judge's decision rate with a chi-square
  test of equal rates, the many-judge slope test of nested judges, and the
  comparison of judges within groups of cases with adjusted p-values and
  slopes; `mmshr()` gives the shares of sentences at, below and above a
  mandatory minimum. `proof("jdgiv")` and `proof("jmono")` name the Lean
  theorems.
* `agemix()` gives the aggregate age-crime curve of a mixture of latent
  types and the covariance across types between ages that a single type
  forces to 0, with its permutation test on panel data; `cumrsk()` gives the
  cumulative risk of a first arrest, conviction or admission by age from a
  life table of first events, with Greenwood's interval or with deaths
  competing; `jpipe()` gives rates between stages of the justice pipeline
  from window sums lagged by the time between the stages. `proof("agemix")`
  names the Lean theorem.
* `mndla()` classifies segregation placements under the Mandela Rules (more
  than 15 days, or Rules 43 and 44 with hours out of cell) and gives the
  qualifying share by group over placements, people with any qualifying
  placement and people by their cumulative days, with exact intervals.
* `stkflw()` relates person-days, the average daily population, admissions
  and the length of stay (any two give the rest) with changes, rates and the
  stock an earlier release removes; `lostay()` summarises lengths of stay,
  their bands and thresholds and the lognormal, exponential and Pareto fits;
  `bklog()` applies Little's law to a docket; `dktflw()` gives clearance,
  rolling clearance, the pending caseload, Little's law on published figures
  and the trend of the outcome mix; `lftab()` builds the actuarial life
  table; `trnmat()` gives the transition matrix between states with its
  stationary distribution. `proof("stkflw")` and `proof("bklog")` name the
  Lean theorems.
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
* `vrfpol()` refuses a CSV with a quote left open in its first rows as a
  file that "could not be read cleanly"; it used to call the file empty.
  A file without a newline at its end no longer lets a warning escape
  from the reading of its header.
* `vrfpol()` ends its refusal of negative exposure values with the word
  it was missing: "a concentration cannot be negative".
* `vrfpol()` refuses a gzip file cut short inside its data as "truncated
  or damaged"; it used to read the part before the cut without a word,
  so the report rested on some of the rows only. The data read are now
  checked against the length and CRC-32 at the end of the file.
* `missgrid()` keeps each letter of the header under its own column when
  the names have different lengths; the letters of a longer name used to
  slide left under the column of a shorter one. The command line had the
  same fault and is fixed with it.
* `sarfc()` refuses an AR(1) fit whose coefficient reaches 1 or beyond with
  "the fitted ARMA part is not stationary.", as it does for higher orders;
  it used to go on with an infinite or negative starting variance.
* A number with no bounds given badly is refused as "must be one number."
  (`beta` of `nomutl()`, `nompr()` and `nomll()`; `c0`, `d0`, `e0`, `a0`
  and `b0` of `dynirt()`) or "must be one number, neither infinite nor
  missing." (`rho` of `dtrid()`, `a` of `inspg()`, `cap` and `floor` of
  `bandx()`); the messages used to name a bound of "-inf" or "-Inf" that
  does not exist. The command line gives the same messages.
* `pedlr()` refuses a pedigree whose tables would pass its size bound
  before it builds them: an untyped child of untyped parents at a locus of
  2000 alleles used to try to allocate 8e9 cells and abort R. The mutation
  matrices are made only with a mutation model, and refused past the same
  bound (more than 8192 alleles); they were made even with `mut = "none"`.
  `gmixlr()` no longer lists the genotypes of unknown contributors when
  there are none, and `ibench()` refuses more than 8190 columns of `x`,
  whose work arrays grow with the square of their number.
* The Bayesian samplers `bamscl()`, `bmdsf()`, `bunfld()`, `cjrirt()`,
  `anomin()` and `ordirt()` refuse an `nsamp` whose draws would pass 2^26
  values before any is kept, and all seven, with `dynirt()`, refuse a
  `burn` plus `nsamp` (times `thin`) past the largest integer; a large
  `nsamp` used to abort R, and the sweep count overflowed. `optcls()`,
  `smacf()`, `nmdsf()`, `mlsmu()`, `smunf()`, `ordocl()`, `indscl()`,
  `dwnom()` and `irtem()` refuse an `ndim` whose configurations would pass
  2^26 values, `anomin()` more dimensions than legislators, `scboot()`
  too many replicates and `anchv()` too many categories, all before
  allocating; `nombt()` no longer keeps its replicates (it draws them
  twice from the seed, for the means and then the deviations), so any
  `nboot` is only long; `smacf()` can now be interrupted. The command line
  gives the same messages.
* `arimfc()` refuses a horizon `h` past 2^26 and `fpdens()` a grid past
  2^26 points before making them (both used to abort R); the forecast
  and the steps of `fpdens()` can now be interrupted. The command line
  gives the same messages.
* `bkrig()` counts the points of its blocks before it lays them out and
  refuses more than 2^26: `nq = .Machine$integer.max` overflowed the count
  and crashed R. The block integrals can now be interrupted within a block.
  The command line gives the same message.
* `evtst()` refuses `leads` plus `lags` past the number of rows, which
  cannot be estimated (each event time needs a row), before it builds the
  indicators; a large window used to abort R. `didfe()` with an `mtch`
  larger than the number of comparison units matches as that many does,
  instead of sizing its tables by `mtch`. The command line does the same.
* `hwksim()` refuses an `nmax` past 2^26 and `sdbsim()` more than 2^26
  `steps`, before the arrays they size are made (both used to abort R).
  `kscan()` compares each simulated maximum as it is drawn instead of
  keeping all `nsim` of them, so a large `nsim` is only long, and can be
  interrupted. The command line gives the same messages.
* `lnroc()` no longer keeps its bootstrap replicates: it draws them twice
  from the seed, summing the areas and then their squared deviations, so
  the errors are those of the kept replicates and a large `nboot` (which
  used to abort R) is only long.
* `inlhist()` refuses more than 2^26 `bins`; the command line gives the
  same message.
* `bayout()` gives the negative binomial bound for large counts (more than
  about 750 a period in the reference window): the search started where
  the probability underflows to 0 and never ended. Counts so large that
  no bound can be found are refused. `evfilt()` stops selecting when no
  eigenvector has a finite criterion (values whose sums overflow); it
  used to add one that does not exist and abort R.
* `dualcont()`, `frmaud()` and the readers of `cnsprof()` results check
  the columns they read even when the class says where the object came
  from, once their other arguments are checked; an object of that class
  with a column missing crashed R.
* `tpdmc()` refuses more than 2^26 `steps` before it makes the path; the
  command line gives the same message.
* `sirex()` checks `level` before it reads it; an empty `level` was read
  past its end and crashed R.
* `tfrp()` refuses a `ymax` of 2^26 or more and `aorwt()` more than 2^26
  bins, before anything is made; `fairsim()` makes as many records as
  asked, or stops with "fairsim: not enough memory for `n` values." The
  command line gives the same messages.
* `crwalk()` makes its simulated paths once, as the matrices it returns
  (they were made twice), and stops with "not enough memory for `h`
  times `nsim` values." when they cannot be made; `vgemp()` refuses more
  than 2^26 distance bins. The command line does the same.
* `lnsize()` refuses a lineup size `k` past 2^26; the command line gives
  the same message.
* `evtst()` refuses more than 2^26 `draws` for the sensitivity interval,
  and `wcbx()` a `reps` whose draws (one value per cluster and outcome)
  would pass 2^26, in R and in the command line. The draws of `evtst()`,
  the bootstrap of `cbatt()`, the placebo replicates of `sdidx()` and the
  permutation and power replicates of `wcbx()` can now be interrupted.
* `enpath()`, `dmlx()` and `suplrn()` refuse an `nlam` past 2^26, in R
  and in the command line.
* `windiv()` refuses a `width` that divides 360 into more than 2^26
  sectors, or into more columns (areas times sectors) than an integer
  counts; a tiny `width` used to overflow the count of columns. The
  command line gives the same messages.
* The robustness iterations of `lwss()`, the replicator steps of
  `insphs()`, the cells of `jtsurf()`, the EM steps of `mcartst()` and
  the bootstrap refits of `dstfit()` can now be interrupted. Loops and
  sizes that multiplied two integer counts are counted in a wide type.
* `exprun()` and `runfx()` refuse a `cap` past 2^26, `jntx()` an `ngrid`
  past 8192 (its grid has `ngrid` squared points) and `dosebin()` an `nq`
  past 2^26. These used to stop with R's own "cannot allocate vector"
  error. The command line gives the same messages where it has the
  function. The cluster bootstrap of `ppmlfe()`, `ivfe()` and `ppmlcf()`
  makes its draws a block at a time instead of holding the clusters by
  `nboot` counts at once, with the same draws in the same order, so any
  number of clusters and draws fits.
* `dmlx()` and `suplrn()` refuse more inner folds (`nfin`) than rows of
  `x`, before their work arrays are sized by it; the command line gives
  the same message.
* `cconc()` refuses an `nsim` past 2^26, `fairdf()` an `nboot` past 2^25
  and `blpal()` more than 2^26 colours, instead of R's "cannot allocate
  vector" error; `cbootmean()` and `cunif()` make results of any length
  R can hold, and stop with "not enough memory for `B` values." (or `n`)
  when it cannot. The command line gives the same messages.
* `stacel()` with a `top` above the cells of its grid keeps them all, as
  that many does, instead of sizing its work arrays by `top`.
* No exception leaves a C entry point of the library: every `clab_*`
  function returns `CLAB_ENOMEM` (new, 4) when an allocation fails in
  it, `clab_jrun()` then answers "fn: not enough memory for this size.",
  and the command-line program prints that and exits with status 1 for a
  failed allocation anywhere, never aborting. Python raises a
  `MemoryError` with the message. The WebAssembly build, which has no
  exceptions, stops the module and the JavaScript binding reports the
  same words (`clab_oomed()`, new) and asks for the module to be loaded
  again. In R the entry points still give the function's name.
* Sizes and products of sizes are 64-bit everywhere they were `long`,
  which is 32 bits on Windows: the work arrays of the learners
  (`enpath()`, `dmlx()`, `suplrn()`, `rfor()`), the neighbour tables, the
  bins of `rddens()` and others. `clab_dmlxsz()` and `clab_stksz()` now
  return `int64_t` and take an `int64_t *`.
* The quantile functions of the command line, the JSON entry and Python
  give R's answers at a point mass: the beta quantile with both shapes 0
  is 1/2 at p = 1/2 (it was 1) and compares p itself in the upper tail as
  R's `qbeta()` does; the gamma and chi-square quantiles with shape or
  degrees of freedom 0 are infinite at p = 1 (they were 0).
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
* Design-based causal effects of a 0/1 treatment. `cbatt()` gives
  Hajek-normalised inverse probability weighting, g-computation,
  subclassification and caliper matching on the propensity (ATE, ATT or
  ATC), and entropy balancing for the ATT, with bootstrap errors;
  `smdbal()` the standardised mean differences and variance ratios before
  and after weighting with the overlap of the propensity; `evlx()` the
  E-value of a risk, odds or hazard ratio or a standardised difference.
* `svymx()` gives design-based means, totals, proportions, ratios and
  variances from stratified, clustered and weighted samples, by domain and
  with differences between domains, with linearisation errors, the
  finite-population correction and a description of the design.
