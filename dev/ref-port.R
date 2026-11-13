# The recorded values of the parity tests (tests/testthat/ref-port.rds):
# the cases of tests/testthat/helper-port.R run against the installed
# rmorie and rmoriebricklayer. From the package root:
#   Rscript dev/ref-port.R [LIB]
# LIB: a library holding rmorie and rmoriebricklayer. The tests use the
# live packages when installed and this file otherwise.
a <- commandArgs(TRUE)
if (length(a)) .libPaths(c(a[1L], .libPaths()))
test_path <- function(...) file.path("tests", "testthat", ...)
sys.source(test_path("helper-port.R"), envir = environment())
ref <- lapply(pcase, function(f) f())
saveRDS(ref, test_path("ref-port.rds"), version = 2)
cat(length(ref), "cases written to", test_path("ref-port.rds"), "\n")
