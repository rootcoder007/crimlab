# The spatial crime methods (near repeats, hot spots, scans, spatial
# regression, point patterns, geographic utilities), whose shapes their own
# test files check; added to spx, which test-crimlab.R's shape test leaves
# out.
spc <- c("knxnr", "kdsurf", "dbscn", "hsper", "ebshr")
spx <- c(spx, spc)
