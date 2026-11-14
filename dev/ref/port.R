# Reference values of the parity tests (tests/testthat/helper-port.R):
# its cases run against rmorie and rmoriebricklayer, written to
# tests/testthat/ref/port.txt. From the package root, with both packages
# in the library:
#   Rscript dev/ref/port.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("rmorie", "rmoriebricklayer")
source(file.path("dev", "ref", "put.R"))
test_path <- function(...) file.path("tests", "testthat", ...)
sys.source(test_path("helper-port.R"), envir = environment())
ref <- lapply(as.list(pcase), function(f) f())
put(ref[sort(names(ref))], "port", need)
