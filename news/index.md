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
