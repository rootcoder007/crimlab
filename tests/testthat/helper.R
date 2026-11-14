# Bad input, table-driven: each case is the message expected, then the
# arguments of one call to f.
#' @srrstats {G5.2b} errs() triggers each refusal and matches its message.
# Every message starts with the name of the function that raised it.
errs <- function(f, cases) {
  for (a in cases) {
    e <- expect_error(do.call(f, a[-1L]), a[[1L]], info = f)
    expect_match(conditionMessage(e), paste0("^", f, ": "), info = f)
  }
}

# Reference values from other packages, stored with the tests: the list
# written to ref/<name>.txt by dev/ref/<name>.R, whose first line names
# the packages and versions that computed it.
refget <- function(name) dget(test_path("ref", paste0(name, ".txt")))
