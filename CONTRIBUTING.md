# Contributing to csls

Bug reports and pull requests are welcome at
<https://github.com/rootcoder007/csls/issues>.

A bug report is most useful with a small reproducible example: the call,
what it returned, and what you expected.

Before a pull request:

* Computation belongs in the C++ core (`src/core.cpp`, declared in
  `src/csls.h`); the R functions check arguments and call it.
* Every change comes with tests in `tests/testthat`. Coverage of R and C++
  must not fall below the floor in `.github/workflows/coverage.yml`.
* `lintr::lint_package()` must report no lints, and
  `R CMD check --as-cran` must pass with no errors or warnings.
* Documentation is written with roxygen2; regenerate `man/` and
  `NAMESPACE` with `roxygen2::roxygenise()`.
* A change to a result that `proof()` names needs the Lean proof in
  `inst/lean` to still build.

By contributing you agree that your contribution is licensed under the
package licence (AGPL-3).
