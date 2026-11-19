# Tests

`tests/testthat.R` runs one test file (`test-crimlab.R`, a few seconds)
unless the environment variable `NOT_CRAN` is `"true"`; with it, it runs
every file under `tests/testthat/`. CRAN therefore runs the short set,
and the project's own checks, which set `NOT_CRAN`, run all of it on
every push, as `devtools::test()` does on a developer's machine.

    NOT_CRAN=true R CMD check crimlab_*.tar.gz

The full suite needs no network, no data beyond the files in
`tests/testthat/`, and no package outside `Suggests`; tests that compare
with another package skip when it is not installed. It takes about five
minutes on one core and runs on Linux, macOS and Windows. The same
functions are also tested from the command line by `cli/test.sh` against
`cli/test.tsv`, the table of cases every language binding shares.
