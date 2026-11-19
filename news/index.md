# Changelog

## crimlab 0.0.2

- The argument that names a method is `mthd` (it was `meth`), and so is
  the part of a result that records it: in
  [`lodimp()`](https://rootcoder007.github.io/crimlab/reference/lodimp.md),
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md),
  [`wcbx()`](https://rootcoder007.github.io/crimlab/reference/wcbx.md),
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md),
  [`cortbl()`](https://rootcoder007.github.io/crimlab/reference/cortbl.md),
  [`topcors()`](https://rootcoder007.github.io/crimlab/reference/topcors.md),
  [`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
  [`scboot()`](https://rootcoder007.github.io/crimlab/reference/scboot.md),
  [`chomap()`](https://rootcoder007.github.io/crimlab/reference/chomap.md),
  [`smmap()`](https://rootcoder007.github.io/crimlab/reference/smmap.md)
  and
  [`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md),
  in R and in every other language. A call that names the old argument
  is refused.

- Plots are black and grey on a white ground with a black frame and a
  bold title, as a journal prints them: the tinted panel and its grid
  are no longer the default. `plotopt(palette = "okabe")` brings the
  colours back, and `plotopt(grid = TRUE, panel = "#EEF1F6")` the panel
  and grid. The `"grey"` palette starts at black and has no grey too
  pale to read on white.

- Intervals are drawn with a cap at each end.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  and
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md)
  with `pal = "diverging"` is blue before adoption and red after.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  power curve of `wcbx(mthd = "power")`.
  `plotevgrid(x, style = "stata")`, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  and
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md)
  with it, draws hollow points joined by a line with dashed intervals.
  The Goodman-Bacon panel of `plot(didfe())` names the two-way estimate
  its line marks.

- New functions for the review of restrictive confinement, after the
  analyses of Sprott and Doob (2021) and of Sprott, Doob and
  Iftene (2021) of Canada’s Structured Intervention Units, in R and in
  every other language.
  [`rvwtim()`](https://rootcoder007.github.io/crimlab/reference/rvwtim.md):
  whether a review due by a set day was held, the stays past that day
  never reviewed, and how long the reviews held took.
  [`entcmp()`](https://rootcoder007.github.io/crimlab/reference/entcmp.md):
  the days of a stay on which an entitlement (hours out of cell, human
  contact) was missed, by how often, with the hours received on those
  days and the part refusals account for.
  [`dmvar()`](https://rootcoder007.github.io/crimlab/reference/dmvar.md):
  whether decision makers decide alike, with each one’s rate, the
  chi-square test, and the test with the small ones and with each one
  left out.

- New plots. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  `wcbx(mthd = "perm")` draws the permuted statistics of each outcome
  with the observed one, as a histogram or a density (the result now
  keeps them, as `null`, in every language).
  `plot(didfe(loo = TRUE), type = "hist")` draws the leave-one-out
  estimates as a histogram with a summary.
  [`plotcoh()`](https://rootcoder007.github.io/crimlab/reference/plotcoh.md)
  draws the mean outcome of each adoption cohort around its adoption,
  and returns the points in every language.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  with `what = "value"` draws the values over the periods.

- The command-line program and the library are built with the arithmetic
  of SSE2 on 32-bit x86, as every other platform computes. With x87
  arithmetic and a distribution’s options (link-time optimisation, no
  frame pointer) three results of the test table differed in their last
  digits, and the openSUSE i586 package of 0.0.1 failed its checks. A
  build that names `-mfpmath` itself keeps its choice.

- A check of CRAN runs one test file, the compiled self-test of every
  kernel and one call of every exported function; the whole suite, which
  takes minutes, runs in the project’s own checks on every push
  (`NOT_CRAN=true`).

- The package installs faster: R compiles its C++ sources in 40 groups
  (a unity build, `src/u_*.cpp`) and no longer one file at a time, so
  the headers are read 40 times, not 337. Nothing a caller sees changes.

- Outside R the plots are drawn too. `plotsvg()` (Python, Go, Rust,
  Java, JavaScript and the command line) runs a plot function on its
  arguments and returns the picture as a standalone SVG document with
  the rows drawn:
  [`plotlrnz()`](https://rootcoder007.github.io/crimlab/reference/plotlrnz.md),
  [`plotgain()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotcdfs()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotrdg()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotvln()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plothex()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md),
  [`plotfnl()`](https://rootcoder007.github.io/crimlab/reference/plotlrnz.md),
  [`plotcal()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plotmult()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plotbump()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plothrz()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md),
  [`plotevgrid()`](https://rootcoder007.github.io/crimlab/reference/plotevgrid.md)
  and
  [`plotcoh()`](https://rootcoder007.github.io/crimlab/reference/plotcoh.md),
  and the results of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md),
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md),
  [`synthc()`](https://rootcoder007.github.io/crimlab/reference/synthc.md)
  and
  [`sdidx()`](https://rootcoder007.github.io/crimlab/reference/sdidx.md),
  with the look of the R plots (black and grey, a black frame, a bold
  title; `palette = "okabe"` for colours, `grid = TRUE` for a grid). One
  drawing code in the C++ core serves every language.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`synthc()`](https://rootcoder007.github.io/crimlab/reference/synthc.md)
  and
  [`sdidx()`](https://rootcoder007.github.io/crimlab/reference/sdidx.md)
  draws the treated path solid and its comparison dashed, as plain
  lines, with a dotted line where the policy was adopted.
  `plot(x, type = "o")` puts a point mark on every period, as before.

- Every package carries the title of the R package, “Statistical Toolkit
  for Criminology and Sociolegal Research”: the Python, JavaScript,
  Rust, Java and Go packages, the command’s manual page, the C header
  and the Homebrew, Scoop, Debian, RPM and container descriptions said
  “criminology and sociolegal statistics”. The README of the Python,
  JavaScript and Rust packages, which is what each registry shows, says
  how to cite crimlab.

- Fixed: `ecolrg(family = "negbin")` on counts with no overdispersion.
  The search for the size sends it off without bound there, and past
  about ten million its score is no more than rounding; on 32-bit x86
  without SSE2 (x87 arithmetic) a step taken from that noise ended at a
  wrong fit, with a lower likelihood. The search now ends at the cap
  once the size passes a million, where the fit is the Poisson fit to
  six digits. Results elsewhere do not change.

- The tests of the command run on x87 as well: its test table compares
  the numbers inside a text result, and two cases whose last digits
  differ there have a bound.

## crimlab 0.0.1

- The results of the Indigenous health and wellbeing functions have
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: gaps
  between groups with their intervals
  ([`indgap()`](https://rootcoder007.github.io/crimlab/reference/indgap.md),
  [`trchlthgap()`](https://rootcoder007.github.io/crimlab/reference/trchlthgap.md)),
  wellness beside distress
  ([`dualcont()`](https://rootcoder007.github.io/crimlab/reference/dualcont.md),
  [`welprof()`](https://rootcoder007.github.io/crimlab/reference/welprof.md),
  [`holbal()`](https://rootcoder007.github.io/crimlab/reference/holbal.md)),
  community indices and Community Well-Being scores
  ([`commidx()`](https://rootcoder007.github.io/crimlab/reference/commidx.md),
  [`cnscwb()`](https://rootcoder007.github.io/crimlab/reference/cnscwb.md),
  [`cwbscr()`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md)),
  standardised prevalence
  ([`stdprev()`](https://rootcoder007.github.io/crimlab/reference/stdprev.md)),
  ordered ratings
  ([`ordcmp()`](https://rootcoder007.github.io/crimlab/reference/ordcmp.md))
  and change between surveys
  ([`healchg()`](https://rootcoder007.github.io/crimlab/reference/healchg.md)).
  A dot plot keeps every row label, in letters sized to the rows, and
  leaves room for its key above the data; lines over periods written as
  text (“2020/21”) are drawn in the order given and labelled as written.

- Open data:
  [`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
  loads a table of the project’s open data service by its key and
  [`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
  lists the tables, with their rows, columns and SHA-256. Everything
  about a table is read from the service’s manifest at the time of the
  call, so a table its publisher replaces is loaded as it is now, with
  no new release. A file is used only when its SHA-256 is the one the
  manifest lists; tables are kept in
  `tools::R_user_dir("crimlab", "cache")` while they still match. `cols`
  picks columns without regard to case or underscores, and `base` (or
  `options(crimlab.data)`) points at a mirror.

- Placement volatility and alert complexity:
  [`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md)
  counts the changes of placement region in each person-year (within a
  record and between consecutive records) and
  [`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md)
  the distinct alert configurations its records span, as defined in
  Ruhela’s study of Ontario’s restrictive confinement records.
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  reduces a release of those records (a data frame, or a key of the open
  data service read by
  [`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md))
  to one row per person and fiscal year, finding the columns without
  regard to case or underscores and taking the fiscal years from the
  data.
  [`cxvol()`](https://rootcoder007.github.io/crimlab/reference/cxvol.md)
  gives, pooled or for each fiscal year, the descriptives, the crude
  rate ratio, the one-to-one propensity matching, the count model with a
  random intercept per regional cluster and the bounds on the effect
  without ignorability;
  [`cxdml()`](https://rootcoder007.github.io/crimlab/reference/cxdml.md)
  the interactive double machine learning estimates (ATE and ATTE);
  [`cxyoy()`](https://rootcoder007.github.io/crimlab/reference/cxyoy.md)
  the year-over-year table of regional changes per 1,000 placements and
  per 1,000 person-years with exact intervals. A single fiscal year with
  a year term, or a sample with no complex or no other person-year, is
  refused in words.
  [`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)
  fits a Poisson or negative binomial count with a normal random
  intercept per cluster by the Laplace approximation, and
  [`cxsim()`](https://rootcoder007.github.io/crimlab/reference/cxsim.md)
  simulates the bias, spread and coverage of the estimators, with or
  without an unobserved confounder.

- [`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
  bounds Cramer’s V right at any total: for a chi-square past about
  1.5e6 (a 2 x 2 table of 1e8 counts) both bounds in R fell below the
  estimate, with warnings from R’s
  [`pchisq()`](https://rdrr.io/r/stats/Chisquare.html), and on the
  command line and in Python the bounds took 149 s at 1e9 counts. R, the
  command line and Python now take them from one routine of the core, an
  integral of the noncentral chi-square whose work does not grow with
  the counts (about a millisecond a bound, to 1e-14 of references of 34
  to 40 digits). With `exact = TRUE` the odds ratio of a 2 x 2 table and
  its interval are summed about their mode: a table of 1e6 counts took a
  minute.

- [`pwranv()`](https://rootcoder007.github.io/crimlab/reference/pwranv.md)
  gives power 1 for a large effect (0.9999999486 at `f = 1000`) and
  returns at once for any `f`; `f = 2^31` and larger never returned.

- [`gstr()`](https://rootcoder007.github.io/crimlab/reference/gstr.md)
  takes the same time whatever the counts: the binomial tail of its
  exposure test summed the terms one by one (20 s for counts near 2e6
  with `nsim = 0`, no end for larger ones).

- [`opteffort()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md)
  gives finite efforts that add to the total for any total: from about
  1e4 the multiplier underflowed and the first cell’s effort was `Inf`.

- [`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md),
  [`olsdx()`](https://rootcoder007.github.io/crimlab/reference/olsdx.md)
  and
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  give the same results for values shifted by a constant (and
  [`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md)
  for values scaled): with values large next to their spread
  [`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md)’s
  z-values moved by 0.4 (1e8) or became `NaN` (1e9),
  [`olsdx()`](https://rootcoder007.github.io/crimlab/reference/olsdx.md)
  moved its slope (1e6) or refused the column as collinear (1e7) and
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  refused it (1e6).

- [`gtwreg()`](https://rootcoder007.github.io/crimlab/reference/gtwreg.md)
  fits an integer `x` with `cons = FALSE` (it stopped in R’s C words).

- [`plothrz()`](https://rootcoder007.github.io/crimlab/reference/plotcal.md)
  refuses a period that is not finite; with one it never returned.

- Refusals and warnings start with the function’s name, in the package’s
  words: about forty functions refused some inputs in base R’s words (a
  date where a column name goes, a function where a vector goes, text
  that is not valid UTF-8, a path in a folder that does not exist, a
  count past what can be allocated, a 1 x 1 matrix for a number), and
  the warnings of
  [`stdprev()`](https://rootcoder007.github.io/crimlab/reference/stdprev.md),
  [`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md),
  [`hetfx()`](https://rootcoder007.github.io/crimlab/reference/hetfx.md)
  and
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md)
  did not name their function. A text argument that is not valid UTF-8
  is refused in R as on the command line (“`ref` holds text that is not
  valid UTF-8”); in Python with the same words.
  [`plotgain()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  no longer says that counts whose total overflows are all 0. A date or
  a time given where the scaling and IRT models take a count (`burn`,
  `nsamp`, `ndim`, `maxit`) is refused; it was read as its days or
  seconds since 1970 and ran for hours.

- [`dktflw()`](https://rootcoder007.github.io/crimlab/reference/dktflw.md)
  takes its clearance intervals from the core’s beta quantile, as the
  command line does: R’s lost its digits, with warnings, once one count
  was 1e13 times the other.

- Proofs: `inst/lean/Volat.lean` (the two measures: bounds, the cases of
  no change and of one record, relabelling, order of the records with a
  counterexample and the sorted case, adding a record, pooling, weighted
  and unweighted means, and two counterexamples on what alert complexity
  counts), `inst/lean/Ident.lean` (identification by adjustment,
  weighting and augmentation, double robustness, exact matching, the
  log-link rate ratio, a negative control, Manski’s bounds and the
  non-identification of the counterfactual arm) and
  `inst/lean/Custody.lean` (Lakner’s average daily population, length of
  stay and expansion estimator with the manual’s worked examples, and
  the arithmetic of year-over-year changes, including what a
  first-to-last column means once a year is added).
  [`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
  lists them.

- [`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
  gives the values the file holds: every column is read as text and
  becomes numbers only when each value is a plain decimal number held
  exactly, so identifiers such as `"007"`, `"7"` and `"0007"` stay
  distinct and a column of `"F"` stays text (it became `FALSE`). A line
  with more or fewer fields than the header is refused, naming the line
  (a long line was wrapped into a row of its own), and so is a table
  whose number of rows is not the one the manifest lists. The manifest
  reader takes all JSON: an empty key, a string of any length (one of a
  few megabytes failed in PCRE), and characters past U+FFFF written as
  surrogate pairs.

- Messages quote invalid UTF-8 bytes of the input as `<xx>`, in every
  refusal and warning of the package and of its C++ binding:
  [`dspbm()`](https://rootcoder007.github.io/crimlab/reference/dspbm.md),
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  and about twenty others carried the bytes themselves, and
  [`substr()`](https://rdrr.io/r/base/substr.html) or
  [`grepl()`](https://rdrr.io/r/base/grep.html) on such a message
  failed.

- The keys that join several columns (the person-years of
  [`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md)
  and
  [`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md),
  the units and days of
  [`expjn()`](https://rootcoder007.github.io/crimlab/reference/expjn.md),
  [`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md),
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  and
  [`updwn()`](https://rootcoder007.github.io/crimlab/reference/updwn.md))
  no longer merge labels holding a carriage return: each value is
  escaped before the values are joined, so no label can imitate the
  join.

- [`qltyrcd()`](https://rootcoder007.github.io/crimlab/reference/qltyrcd.md)
  reads its flags once: a vector of a million flags took minutes,
  copying a row per value, and is now read in a fraction of a second.

- Arguments of the wrong type are refused in the package’s words, not
  base R’s:
  [`gtwreg()`](https://rootcoder007.github.io/crimlab/reference/gtwreg.md)
  takes `x = NULL` (a model with the constant only, which stopped on its
  column names) and refuses a function or an environment as `x`;
  [`grdbin()`](https://rootcoder007.github.io/crimlab/reference/grdbin.md)
  a `col` that is not one name;
  [`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md)
  an `imported` that is not a vector;
  [`dagcol()`](https://rootcoder007.github.io/crimlab/reference/dagcol.md)
  a `triple` that is not labels;
  [`lrcce()`](https://rootcoder007.github.io/crimlab/reference/lrcce.md)
  an `x` that names no column; the plots a `col` that is not colours (a
  function was warned about, a data frame failed in the graphics).

- The help of
  [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  says what it refuses with an R error: only a `poll` that is not one
  pollutant name, as the command line refuses it; it said the function
  never stops with an R error.

- [`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md)
  gives R’s [`qnbinom()`](https://rdrr.io/r/stats/NegBinomial.html)
  bound for every reference sum, in R, the command line and Python. For
  sums near 4,800 (a week of about 800) the probability of 0 was a
  subnormal number: the bound came out too low (805 for 856, raising
  alarms wrongly) and for nearby sums the search never ended
  (`bayout(1:810)` did not return). The quantile is now a search of at
  most about 160 evaluations of the distribution function, and the
  Poisson quantile of the command line and Python takes the same search.

- [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  is in `clab_jrun()`, the command line and Python, with the report R
  gives: the demo data are R’s draws, and the exposure file is read as R
  reads plain and gzip UTF-8 text. A bzip2 or xz file and UTF-16, UTF-32
  or Latin-1 text, which R decodes, are an error report outside R.

- New:
  [`tabrpt()`](https://rootcoder007.github.io/crimlab/reference/tabrpt.md),
  one call on a published table:
  [`catchk()`](https://rootcoder007.github.io/crimlab/reference/catchk.md)
  on its groups,
  [`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
  on the change between periods,
  [`grprr()`](https://rootcoder007.github.io/crimlab/reference/grprr.md)
  on the rates and
  [`trndrr()`](https://rootcoder007.github.io/crimlab/reference/trndrr.md)
  on each group’s trend.

- New:
  [`ovrchk()`](https://rootcoder007.github.io/crimlab/reference/ovrchk.md)
  flags the cells of a coded table of oversight reports that are not of
  the form of their field (a real date, a count, a flag) or of a pattern
  such as a case number.

- New:
  [`ovrkwd()`](https://rootcoder007.github.io/crimlab/reference/ovrkwd.md)
  marks the reports that mention the words of keyword groups, whole
  words in any case, with a preset for mental health and race and any
  groups of one’s own.

- New:
  [`ovrsoc()`](https://rootcoder007.github.io/crimlab/reference/ovrsoc.md)
  counts the subject officials of an oversight report by the first of
  six rules that applies (a roster, numbered tags, a plural, a singular
  mention, a statement of none) and says which rule and words gave the
  count.

- New:
  [`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
  codes police oversight reports into the fields of a schema
  ([`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md)),
  by rules that are deterministic and offline, one row per report.

- New:
  [`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md),
  the field schemas of oversight reports: each field read by a rule
  (text after a cue, a number, a date in the cue’s sentence, numbered
  tags, a flag and others), with a preset for the director’s reports of
  Ontario’s Special Investigations Unit.

- New:
  [`ovrdte()`](https://rootcoder007.github.io/crimlab/reference/ovrdte.md)
  reads the first date written out in a text, in English or French, as
  `YYYY-MM-DD`.

- New:
  [`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md),
  the plain text of a police oversight report’s page, with each block on
  a line of its own so that section headings stand alone; the first of
  the tools that code oversight reports.

- New:
  [`synfix()`](https://rootcoder007.github.io/crimlab/reference/synfix.md)
  draws example data in the shape of a published table: presets with the
  columns of Ontario’s corrections tables (segregation, deaths in
  custody, use of force), or any table given as the values of its
  columns. The draws are R’s after
  [`set.seed()`](https://rdrr.io/r/base/Random.html), in every language.

- New:
  [`pscred()`](https://rootcoder007.github.io/crimlab/reference/pscred.md),
  the credit for custody before sentence under s. 719(3) and (3.1) of
  the Criminal Code, at up to one and a half days a day, and the time
  left to serve.

- New:
  [`olsdx()`](https://rootcoder007.github.io/crimlab/reference/olsdx.md),
  least squares with the checks Algahtany, Kumar and Barclay (2022)
  report for their weather and crime regressions (standardised
  coefficients, tolerance and variance inflation, the Durbin-Watson
  statistic, leverage, studentised residuals and Cook’s distance) and
  the expected count, or a new month’s, at given weather.

- Preset `laniyonu-gent-cnt`: the tract layout of Laniyonu (2018), stops
  on the three-level gentrification factor (the `gent` recipe) with the
  population as the offset and the crime and 311 controls, for
  [`cntreg()`](https://rootcoder007.github.io/crimlab/reference/cntreg.md).

- New:
  [`dagcol()`](https://rootcoder007.github.io/crimlab/reference/dagcol.md),
  the colliders of a directed acyclic graph given by its edges, shielded
  or not (v-structures), and whether a given triple is one.

- New:
  [`tpsurf()`](https://rootcoder007.github.io/crimlab/reference/tpsurf.md),
  a thin-plate spline surface through values at points in the plane
  (Duchon 1977), smoothed by `lam`, on a grid with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

- New:
  [`eofst()`](https://rootcoder007.github.io/crimlab/reference/eofst.md),
  the empirical orthogonal functions of a matrix of periods by places:
  spatial patterns, their series over time and the share of the variance
  of each.

- New:
  [`prpmap()`](https://rootcoder007.github.io/crimlab/reference/prpmap.md),
  a proportional symbol map: circles whose areas are proportional to the
  values, drawn largest first.

- New:
  [`ptmap()`](https://rootcoder007.github.io/crimlab/reference/ptmap.md),
  a point map in base graphics, with clusters from
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md)
  or given, each in its own colour and shape with its number at its
  centre;
  [`dbscn()`](https://rootcoder007.github.io/crimlab/reference/dbscn.md)
  results now have a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
  draws them through it.

- New:
  [`smmap()`](https://rootcoder007.github.io/crimlab/reference/smmap.md),
  small multiples of a choropleth map, one panel per period with the
  breaks and colours shared by every period.

- New:
  [`chomap()`](https://rootcoder007.github.io/crimlab/reference/chomap.md),
  a choropleth map in base graphics from polygon coordinates, with no
  spatial package: quantile, equal or given class breaks computed in the
  compiled kernel.

- New:
  [`gocoef()`](https://rootcoder007.github.io/crimlab/reference/gocoef.md),
  the slopes of a
  [`gologt()`](https://rootcoder007.github.io/crimlab/reference/gologt.md)
  fit laid out by threshold, with their standard errors and Wald tests.

- New:
  [`gologt()`](https://rootcoder007.github.io/crimlab/reference/gologt.md),
  the generalised ordered logit (Peterson and Harrell 1990; Williams
  2006): an intercept and slopes for each threshold of an ordered
  outcome, chosen covariates held parallel, with the Wald test of each
  covariate’s slopes and the likelihood-ratio test against proportional
  odds.

- New:
  [`ulcate()`](https://rootcoder007.github.io/crimlab/reference/ulcate.md),
  the U-learner of conditional average treatment effects (Kunzel et
  al. 2019) on the cross-fitted learners of
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
  with a final learner of its own and a clipped denominator.

- New:
  [`xpshap()`](https://rootcoder007.github.io/crimlab/reference/xpshap.md),
  the sampling Shapley values of Strumbelj and Kononenko (2014) for
  single predictions, with their standard errors.

- New:
  [`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md),
  the accumulated local effects of one input (Apley and Zhu 2020), with
  the binning and centring of their ‘ALEPlot’.

- New:
  [`xpcp()`](https://rootcoder007.github.io/crimlab/reference/xpcp.md),
  ceteris paribus profiles: the prediction for single rows as one input
  moves over a grid.

- New:
  [`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md),
  the partial dependence of a prediction on one input (Friedman 2001)
  over an even or a quantile grid.

- New:
  [`xpimp()`](https://rootcoder007.github.io/crimlab/reference/xpimp.md),
  the permutation importance of a risk tool’s inputs (Breiman 2001;
  Fisher, Rudin and Dominici 2019), for an R function or the learners of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md),
  with the shuffles drawn in the compiled core and returned.

- New:
  [`mipmm()`](https://rootcoder007.github.io/crimlab/reference/mipmm.md),
  multiple imputation of numeric columns by predictive mean matching in
  chained equations, with an analysis pooled by Rubin’s rules and the
  degrees of freedom of Barnard and Rubin, after
  [`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md).

- New:
  [`rscchi()`](https://rootcoder007.github.io/crimlab/reference/rscchi.md),
  the Rao-Scott test of independence in a two-way table from a
  stratified, clustered and weighted survey, first and second order and
  as an F, with the generalised design effects, beside
  [`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
  and
  [`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md).

- New:
  [`rcttml()`](https://rootcoder007.github.io/crimlab/reference/rcttml.md),
  targeted maximum likelihood for a randomised trial with covariate
  adjustment (Moore and van der Laan 2009), for continuous and binary
  outcomes, beside the g-computation and the difference in means.

- New:
  [`trnsp()`](https://rootcoder007.github.io/crimlab/reference/trnsp.md),
  a trial’s average effect transported to a target population, or
  generalised to the population it came from, by the outcome model and
  by inverse odds weighting, with bootstrap errors.

- New:
  [`gensw()`](https://rootcoder007.github.io/crimlab/reference/gensw.md),
  the sampling score of a trial against its target population and the
  weights that generalise or transport its units, with the logit
  standardised difference, Tipton’s generalisability index, the
  effective size and the covariate balance.

- New:
  [`specrv()`](https://rootcoder007.github.io/crimlab/reference/specrv.md),
  the specification curve of Simonsohn, Simmons and Nelson (2020): the
  treatment’s estimate for every outcome and every subset of optional
  covariates, with the joint inference of the median, the count
  significant and the summed t against permuted or bootstrapped null
  curves.

- New:
  [`plntfx()`](https://rootcoder007.github.io/crimlab/reference/plntfx.md),
  the recovery of planted effects: the outcome shuffled, an effect of
  known size added to the treated and the analysis run again, reporting
  bias, spread, coverage and rejection.

- New:
  [`falsctl()`](https://rootcoder007.github.io/crimlab/reference/falsctl.md),
  falsification with negative-control outcomes: the treatment’s
  coefficient on outcomes it cannot affect, with classical or HC1
  errors, Kennedy’s permutation p-values and an adjustment for the
  number of controls.

- New:
  [`pwrsim()`](https://rootcoder007.github.io/crimlab/reference/pwrsim.md),
  the power of a two-arm trial by simulation, for normal, binary or
  count outcomes, unequal arms and Welch’s or the permutation test, with
  its Monte Carlo error.

- New:
  [`pwranv()`](https://rootcoder007.github.io/crimlab/reference/pwranv.md),
  the power of the one-way analysis of variance from the noncentral F
  distribution, or the group size that reaches a given power.

- New:
  [`blkprm()`](https://rootcoder007.github.io/crimlab/reference/blkprm.md),
  the permutation test of a treatment effect with the assignment
  permuted within blocks, by enumeration or by draws, for the difference
  in means or the block-weighted difference.

- New:
  [`rndalc()`](https://rootcoder007.github.io/crimlab/reference/rndalc.md),
  randomisation lists by simple, permuted-block or stratified block
  randomisation from a seed, the same list in R and on the command line.

- [`plotrdg()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  and
  [`plotvln()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  take their densities from the compiled kernel that the command line
  and the other languages use, so the curve is the same everywhere and
  on every version of R. The kernel agrees with
  [`stats::density()`](https://rdrr.io/r/stats/density.html) of R 4.3
  and later to rounding.

- New:
  [`cemtch()`](https://rootcoder007.github.io/crimlab/reference/cemtch.md),
  coarsened exact matching (Iacus, King and Porro 2012) with the weights
  for the effect on the treated and the L1 imbalance, as Sampson and
  Winter (2018) match lead-exposed children.

- New:
  [`lrcce()`](https://rootcoder007.github.io/crimlab/reference/lrcce.md),
  the CS-ARDL long-run effects and adjustment speed from the unit
  estimates of
  [`ccepnl()`](https://rootcoder007.github.io/crimlab/reference/ccepnl.md),
  as Ditzen (2021) sets them out.

- New:
  [`ifesc()`](https://rootcoder007.github.io/crimlab/reference/ifesc.md),
  the generalized synthetic control of Xu (2017) on the interactive
  fixed effects model of Bai (2009), with staggered adoption, the number
  of factors by cross-validation and a bootstrap of units.

- New: `ppmlcf(family = "nb")`, the control function with a negative
  binomial second stage and a cluster bootstrap of both stages, as
  Feigenbaum and Muller (2016) estimate lead’s effect on homicide.

- New: `ivfe(se = "dk")`, the standard errors of Driscoll and Kraay
  (1998), robust to correlation across units and over time, as Batkeyev
  and DeRemer (2023) report them.

- New:
  [`idwint()`](https://rootcoder007.github.io/crimlab/reference/idwint.md)
  carries monitor readings to places by inverse-distance weighting of
  the nearest monitors and averages them by zone, as Rahman et
  al. (2023), Wongchanapai et al. (2025) and Wesselbaum (2022) assign
  exposure;
  [`gtwreg()`](https://rootcoder007.github.io/crimlab/reference/gtwreg.md)
  is the geographically and temporally weighted regression of Kuo and
  Putra (2021), its bandwidth and time scale chosen by cross-validation
  or AICc;
  [`asrule()`](https://rootcoder007.github.io/crimlab/reference/asrule.md)
  gives the association rules between two binned variables that
  Wongchanapai et al. (2025) mine with Apriori.

- The command line computes the quantiles and medians of
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md),
  `ccross(ns =)` and
  [`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md)
  as R’s [`quantile()`](https://rdrr.io/r/stats/quantile.html) and
  [`median()`](https://rdrr.io/r/stats/median.html) do, and the
  restricted cubic spline’s cube as a product in both, so that a spline
  surface of
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md)
  is the same fit in R and outside it (it differed by about 1e-8).

- Fixed, after a review of the checks:
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md)
  and
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  with `se = "twoway"` and other than two fixed effects,
  [`rddens()`](https://rootcoder007.github.io/crimlab/reference/rddens.md)
  and
  [`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)
  with fewer than 10 rows,
  [`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)
  with a constant treatment, a refused row of
  [`polbrdarea()`](https://rootcoder007.github.io/crimlab/reference/polbrdarea.md)
  and an outcome
  [`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
  lacks for PM2.5 are refused in words that name what is wrong;
  [`tailwd()`](https://rootcoder007.github.io/crimlab/reference/tailwd.md)
  refuses a negative speed and `dwnexp(lonlat = TRUE)` a latitude past
  90;
  [`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md)
  gives winds that cancel, and calm, no direction (`NA`); a whole number
  past 2147483647 is refused naming that limit; and the bin labels of
  [`invsn()`](https://rootcoder007.github.io/crimlab/reference/invsn.md),
  [`dosebin()`](https://rootcoder007.github.io/crimlab/reference/dosebin.md)
  and
  [`updwn()`](https://rootcoder007.github.io/crimlab/reference/updwn.md)
  write `Inf` without a space, as the command line does.

- Fixed: Python
  [`orrisk()`](https://rootcoder007.github.io/crimlab/reference/orrisk.md)
  refused `or_`, although an argument named as a Python keyword takes a
  final underscore elsewhere (`from_`, `lambda_`).

- Fixed:
  [`inlhist()`](https://rootcoder007.github.io/crimlab/reference/inlhist.md)
  lost a value when the range of the values overflowed
  (`c(-1e308, 1e308)`); the breaks are now built from halves there.

- [`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md)
  with very narrow sectors builds its column names in less than half the
  time: each edge is formatted once.

- [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  takes a day column of ISO date text or of days since 1970-01-01, as
  the command line gives dates, and
  [`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md)
  reads ISO date text as dates, in R and on the command line alike;
  before, R’s
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  refused numbers and the command line could not give
  [`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md)
  dates.

- Fixed:
  [`irtem()`](https://rootcoder007.github.io/crimlab/reference/irtem.md),
  [`optcls()`](https://rootcoder007.github.io/crimlab/reference/optcls.md),
  [`dwnom()`](https://rootcoder007.github.io/crimlab/reference/dwnom.md),
  [`cjrirt()`](https://rootcoder007.github.io/crimlab/reference/cjrirt.md),
  [`nomll()`](https://rootcoder007.github.io/crimlab/reference/nomll.md),
  [`irtlik()`](https://rootcoder007.github.io/crimlab/reference/irtlik.md)
  and
  [`nombt()`](https://rootcoder007.github.io/crimlab/reference/nombt.md)
  accepted votes other than 1, 0 and `NA`, and a matrix of no votes;
  they refuse them now in the words of
  [`anomin()`](https://rootcoder007.github.io/crimlab/reference/anomin.md).
  [`ordocl()`](https://rootcoder007.github.io/crimlab/reference/ordocl.md)
  refuses a matrix of no responses, and a vector at once.

- Fixed:
  [`clncolnms()`](https://rootcoder007.github.io/crimlab/reference/clncolnms.md)
  could return repeated names (`c("x", "x", "x_2")` gave `"x_2"` twice);
  a number whose name is taken is now skipped. A tab inside a name is a
  word break, as the help says, and no longer deleted.

- Fixed:
  [`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md)
  refused columns in different units (income next to a rate) as
  collinear; the rank check and the EM step’s eigenvalue floor are now
  on the scale of the correlations.

- Fixed: `rdloc(bw = "ik")` and
  [`rddens()`](https://rootcoder007.github.io/crimlab/reference/rddens.md)
  refused a running variable in metres or in thousands of kilometres;
  the bandwidths now scale with the running variable, and the estimates
  do not change.

- Fixed:
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  and the other least-squares fits refused columns about 1e7 apart in
  scale as collinear; a column is now judged against its own size, and a
  column the fixed effects absorb is found by an explicit check in
  [`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md),
  [`wcbx()`](https://rootcoder007.github.io/crimlab/reference/wcbx.md)
  and the tests of
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md).

- Fixed:
  [`plot.lnkpr()`](https://rootcoder007.github.io/crimlab/reference/plotm.md)
  returned an ROC curve of `NA` and an area of `NA` when some crimes
  were unsolved; the pairs with an unsolved crime are now left out, as
  [`lnkpr()`](https://rootcoder007.github.io/crimlab/reference/lnkpr.md)
  leaves them out of its fit, and the area is the model’s AUC.

- [`polrpttxt()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  is in `clab_jrun()`, the command line and Python: it writes the text
  of a
  [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  report given as JSON, its numbers as R’s
  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) and
  [`format()`](https://rdrr.io/r/base/format.html) write them.
  [`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
  itself runs in R only.

- [`yoypdf()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md)
  is in `clab_jrun()`, the command line and Python, writing the PDF R
  writes; without `file` the PDF is returned as text, and
  [`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md)
  there writes PDF too.

- [`yoyhtml()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md)
  is in `clab_jrun()`, the command line and Python, writing the page R
  writes; without `file` the page is returned as text, and
  [`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md)
  there writes HTML too.

- [`prechk()`](https://rootcoder007.github.io/crimlab/reference/prechk.md)
  is in `clab_jrun()`, the command line and Python, with the publishers’
  codes of `codes.tsv` built into the command line as the presets are
  (`cli/codes.sh` writes `cli/codes.inc`, and CI checks it is up to
  date).

- [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md)
  is in `clab_jrun()`, the command line and Python, where no session
  keeps options: it checks the options given as R does and returns them
  with the defaults for the others.

- [`plotevgrid()`](https://rootcoder007.github.io/crimlab/reference/plotevgrid.md)
  is in `clab_jrun()`, the command line and Python: its `x` maps each
  panel to a result of
  [`evtst()`](https://rootcoder007.github.io/crimlab/reference/evtst.md)
  or
  [`stgdid()`](https://rootcoder007.github.io/crimlab/reference/stgdid.md),
  or to the rows `evtm`, `est`, `lo` and `hi` themselves.

- [`plothex()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  returns its hexagons in order of row, then column, in every locale (it
  ordered them by [`table()`](https://rdrr.io/r/base/table.html) of text
  keys, which sorts in the session’s collation), and is in
  `clab_jrun()`, the command line and Python with the same order.

- [`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md)
  is in `clab_jrun()`, the command line and Python as the conversion it
  makes: the caller gives the wall and CPU seconds (and, where measured,
  the machine’s CPU accounting and the RAPL counters before and after)
  and gets R’s result without `value`. Timing an expression stays in R.

- [`hwpwr()`](https://rootcoder007.github.io/crimlab/reference/hwpwr.md)
  is in `clab_jrun()`, the command line and Python. The result holds the
  rows R returns, the interval `dt` and the last reading as text
  (`snap`), which a later call takes as `since`; R keeps the last two as
  attributes.

- [`plotrdg()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  and
  [`plotvln()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  are in `clab_jrun()`, the command line and Python, returning the rows
  R returns (without the quartiles
  [`plotvln()`](https://rootcoder007.github.io/crimlab/reference/plotrdg.md)
  keeps as an attribute). Their densities come from a kernel that
  reproduces [`stats::density()`](https://rdrr.io/r/stats/density.html):
  the bandwidth rules nrd0, nrd, ucv, bcv and SJ, the binning, the FFT
  convolution and the interpolation, equal to R’s to rounding. In R both
  still call [`stats::density()`](https://rdrr.io/r/stats/density.html).

- New:
  [`dwnexp()`](https://rootcoder007.github.io/crimlab/reference/dwnexp.md)
  builds downwind-exposure instruments from sources of pollution and the
  wind: the sector and transport-wedge rules of Bao and Da (2026) for
  smoke plumes that move from day to day, and the angle and distance
  instrument of Ayesh (2023) for fixed sources.
  [`tailwd()`](https://rootcoder007.github.io/crimlab/reference/tailwd.md)
  gives the tailwind component along a flight course, from which
  Hener (2022) builds the instrument of the landing direction for
  aircraft noise.

- New:
  [`invsn()`](https://rootcoder007.github.io/crimlab/reference/invsn.md)
  takes `by` and `rule` to turn readings into one inversion per day, as
  Batkeyev and DeRemer (2023) require every three-hourly reading of a
  day to be inverted.

- The help pages of
  [`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md),
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
  [`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md),
  [`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md)
  and
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  cite the designs of Karamik and von Graevenitz (2022), Sarmiento
  (2020), Batkeyev and DeRemer (2023), Hener (2022), Bao and Da

  2026. and Rahman et al. (2023), with the arguments that reproduce
        each.

- New:
  [`dirsim()`](https://rootcoder007.github.io/crimlab/reference/dirsim.md)
  draws compositions from a Dirichlet distribution (a gamma variate per
  part over its row sum), from R’s generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) and
  [`rgamma()`](https://rdrr.io/r/stats/GammaDist.html) give the same
  draws; for testing a compositional analysis such as
  [`clrilr()`](https://rootcoder007.github.io/crimlab/reference/clrilr.md)
  before real measurements exist.

- New:
  [`evlx()`](https://rootcoder007.github.io/crimlab/reference/evlx.md)
  takes `true`, a true value other than the null to set the estimate
  against, as the argument of that name in the ‘EValue’ package.

- Fixed: the second sweep found 25 more functions
  ([`lodimp()`](https://rootcoder007.github.io/crimlab/reference/lodimp.md),
  [`lnroc()`](https://rootcoder007.github.io/crimlab/reference/lnroc.md),
  [`ibench()`](https://rootcoder007.github.io/crimlab/reference/ibench.md),
  [`synexp()`](https://rootcoder007.github.io/crimlab/reference/synexp.md),
  [`bndxw()`](https://rootcoder007.github.io/crimlab/reference/bndxw.md),
  [`rptvic()`](https://rootcoder007.github.io/crimlab/reference/rptvic.md),
  [`etohlr()`](https://rootcoder007.github.io/crimlab/reference/etohlr.md),
  [`mndla()`](https://rootcoder007.github.io/crimlab/reference/mndla.md),
  [`crchyp()`](https://rootcoder007.github.io/crimlab/reference/crchyp.md),
  [`clrilr()`](https://rootcoder007.github.io/crimlab/reference/clrilr.md),
  …) whose compiled readers asked for the length of `NULL`, a function
  or an environment; they now refuse it in words.

- Fixed, after a second sweep: a function given as `est` to the exposure
  designs, `end` of a point
  [`datumpmap()`](https://rootcoder007.github.io/crimlab/reference/ltrlrng.md),
  a list for the data of
  [`expconcidx()`](https://rootcoder007.github.io/crimlab/reference/expconcidx.md),
  [`exprspplr()`](https://rootcoder007.github.io/crimlab/reference/exprspplr.md)
  and
  [`polbrdarea()`](https://rootcoder007.github.io/crimlab/reference/polbrdarea.md),
  `ref` and a single `p2` of
  [`jntx()`](https://rootcoder007.github.io/crimlab/reference/jntx.md),
  a bootstrap `seed` past an integer in the exposure estimators, an `x`
  of
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  that is not numbers, an
  [`onlsmry()`](https://rootcoder007.github.io/crimlab/reference/onlsmry.md)
  accumulator without its state, an infinite `cnts` of
  [`bandx()`](https://rootcoder007.github.io/crimlab/reference/bandx.md),
  a non-numeric `source` of
  [`advdif2d()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
  an object without a `crs` in
  [`crsof()`](https://rootcoder007.github.io/crimlab/reference/crsof.md),
  invalid text in
  [`pbbnd()`](https://rootcoder007.github.io/crimlab/reference/pbbnd.md)
  and a table of no rows in
  [`yoysmry()`](https://rootcoder007.github.io/crimlab/reference/yoysmry.md)
  used to fail in R’s words; each is now refused in the package’s, and
  [`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)
  names itself where it named
  [`plcday()`](https://rootcoder007.github.io/crimlab/reference/plcday.md).

- Fixed:
  [`kscan()`](https://rootcoder007.github.io/crimlab/reference/kscan.md)’s
  Monte Carlo p-values used to depend on the order of the rows for a
  given seed (the times were shuffled over the rows as given); the
  shuffles now run over the incidents in the order of their coordinates
  and times. The p-values for a given seed change once; the clusters do
  not.

- Fixed:
  [`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)
  of a table of zeros used to print “n = NA” and give no statistics;
  such a table has no total, and it is now refused (here and on the
  command line).

- Fixed, one function at a time:
  [`bufas()`](https://rootcoder007.github.io/crimlab/reference/bufas.md)
  multiplied a factor or text `r` for its default `buf` before checking
  it (a warning or “non-numeric argument to binary operator”);
  [`gtmp()`](https://rootcoder007.github.io/crimlab/reference/gtmp.md)
  tested a function `q` for `NA` with a warning; the `rgnmap*()` checks
  turned a function into text;
  [`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
  [`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)
  and
  [`expgmm()`](https://rootcoder007.github.io/crimlab/reference/expgmm.md)
  converted a `maxit` past an integer with a warning;
  [`yoypdf()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md)
  did not check `lag`;
  [`dstnskch()`](https://rootcoder007.github.io/crimlab/reference/dstnskch.md),
  [`rsvsmp()`](https://rootcoder007.github.io/crimlab/reference/rsvidx.md),
  [`relbl()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
  [`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md)
  and
  [`grdlvls()`](https://rootcoder007.github.io/crimlab/reference/grdlvls.md)
  took a function for their values or labels;
  [`crfno2()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md),
  [`crfpm25()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md)
  and
  [`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
  failed quoting a function as an outcome or pollutant. Each now refuses
  in its own words.

- Fixed: a colour index past the range of an integer (`col = 1e300`,
  `Inf`) used to pass the plots’ colour check with a coercion warning,
  or fail in [`col2rgb()`](https://rdrr.io/r/grDevices/col2rgb.html)
  (“numerical color values must be \>= 0”), in 12 plot functions; it is
  now refused as any bad colour (“`col` must be colours.”). A 1 x 1
  matrix for a number such as `bins` used to be accepted with a
  recycling warning; it is now refused.

- Fixed: a file path longer than a path can be (1e6 characters) used to
  reach [`path.expand()`](https://rdrr.io/r/base/path.expand.html) and
  its warning (“expanded path length … would be too long”) in
  [`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md),
  [`raplok()`](https://rootcoder007.github.io/crimlab/reference/raplok.md),
  [`wrtxtfb()`](https://rootcoder007.github.io/crimlab/reference/wrtxtfb.md)
  and the `yoy*()` writers; it is now refused, at most 4096 bytes.

- Fixed: text that is not valid UTF-8 given for an argument used to stop
  49 functions inside base R’s string functions (“invalid multibyte
  string”, “input string 1 is invalid UTF-8”) or, in a title, at the
  graphics device; the shared text check now refuses it first, naming
  the argument
  ([`bandx()`](https://rootcoder007.github.io/crimlab/reference/bandx.md),
  [`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md),
  [`gplume()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
  [`fylbl()`](https://rootcoder007.github.io/crimlab/reference/yoylbl.md),
  the `yoy*()` writers, the descriptive statistics, the plots’ `main`,
  …), and a refusal that quotes such text shows its bytes.

- Fixed: `NULL`, a function or an environment given for a vector
  argument used to reach the compiled code and fail there in R’s words
  (“LENGTH or similar applied to NULL object”, “cannot coerce type
  ‘builtin’ to vector”), in 31 functions
  ([`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md),
  [`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md),
  [`tmpaud()`](https://rootcoder007.github.io/crimlab/reference/tmpaud.md),
  [`lodimp()`](https://rootcoder007.github.io/crimlab/reference/lodimp.md),
  [`trnmat()`](https://rootcoder007.github.io/crimlab/reference/trnmat.md),
  [`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md),
  …); it is now refused in the package’s (“`area` must be a vector, got
  NULL.”).

- New:
  [`ccepnl()`](https://rootcoder007.github.io/crimlab/reference/ccepnl.md)
  estimates dynamic heterogeneous panels as Chudik and Pesaran (2015)
  set out: `ylag` lags of the outcome among the regressors, `lags` lags
  of the cross-section averages (by default the integer part of
  `T^(1/3)` when `ylag > 0`), the mean group estimator and its variance,
  and the half-panel jackknife (`bias = "jack"`) or recursive mean
  adjustment (`bias = "rma"`), with the uncorrected estimate alongside
  in `mg0`. The pooled estimator is not made for dynamic models
  (`pool`), as their theory is that of the mean group estimator. Calls
  without the new arguments give the results they gave.

- New:
  [`ccepnl()`](https://rootcoder007.github.io/crimlab/reference/ccepnl.md)
  reports the heterogeneity-bias test of the pooled estimator of Juodis
  and Reese (2026, `hb`) and the number of factors the averages identify
  by their two regularized procedures (`nfac`), and offers the fixed-T
  variance of the pooled estimator of Brown, Schmidt and Wooldridge
  (2026, `pvar = "bsw"`). The help page sets out the five lessons of
  Juodis and Reese and what each asks of a user.

- New:
  [`gmmpnl()`](https://rootcoder007.github.io/crimlab/reference/gmmpnl.md)
  returns the residuals in levels (`resid`), for Pesaran’s test of
  cross-sectional dependence and Moran’s I on each year, and the joint
  test as xtabond2 prints it (`xwald`); its help page maps each estimate
  and test of Andreatta et al. (2023) to a call and has an example of
  their specification on simulated data.

- Changed: `gmmpnl(small = TRUE)` now makes the small-sample corrections
  as xtabond2 3.7.2 computes them (Roodman 2009): the conventional
  one-step variance and `sigma2` times `w / (w - k)`, the robust and
  two-step variances times `(n - 1) / (n - k) G / (G - 1)`, and t and F
  on `n - q - c` or `G - c` degrees of freedom; it used `n / (n - k)`
  throughout. The joint test counts its degrees of freedom as the
  program does, with or without `small`. `small = TRUE` is refused with
  one group, or with no more transformed rows than coefficients.

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
