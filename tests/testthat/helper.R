# Bad input, table-driven: each case is the message expected, then the
# arguments of one call to f.
errs <- function(f, cases) {
  for (a in cases) expect_error(do.call(f, a[-1L]), a[[1L]], info = f)
}
