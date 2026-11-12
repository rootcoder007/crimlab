# Changelog

## crimlab 0.0.0.9000

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
