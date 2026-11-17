# Shared by the placement volatility tests (test-plcvol.R, test-otpyr.R,
# test-cntmix.R, test-cxvol.R, test-cxsim.R).

# print(x) writes a header matching pat and gives x back invisibly.
prt <- function(x, pat) {
  v <- NULL
  expect_output(v <- withVisible(print(x)), pat)
  expect_false(v$visible)
  expect_identical(v$value, x)
}
