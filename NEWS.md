# crimlab 0.0.0.9000

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
