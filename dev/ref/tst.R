# Runs the tests that read reference values. With CLAB_REF set they read the
# regenerated values in that folder instead of tests/testthat/ref/. Run from
# the package root with crimlab installed: Rscript dev/ref/tst.R
f <- paste0("^(spw|sp|spec|kdsurf|cntreg|didfe|ctabx|net|corr|disp|pol|",
            "fairdf|rates|clprep|ref)$")
r <- as.data.frame(testthat::test_dir("tests/testthat", filter = f,
  package = "crimlab", load_package = "installed", reporter = "summary",
  stop_on_failure = FALSE))
cat("pass fail skip error =", sum(r$passed), sum(r$failed), sum(r$skipped),
    sum(r$error), "\n")
if (sum(r$failed) + sum(r$error) > 0) quit(status = 1)
