# Bad input, table-driven: each case is the message expected, then the
# arguments of one call to f.
# Every message starts with the name of the function that raised it.
errs <- function(f, cases) {
  for (a in cases) {
    e <- expect_error(do.call(f, a[-1L]), a[[1L]], info = f)
    expect_match(conditionMessage(e), paste0("^", f, ": "), info = f)
  }
}
