# The exported spatial (spf) and dossier (dsx) functions, whose shapes
# test-sp.R and test-doss.R check; spx is both, left out of the shape test
# of test-crimlab.R.
spf <- c("hvdist", "eqproj", "eqinv", "cgstat", "cgmedn", "webpt", "sdell",
         "cvhull", "nnidx", "ripk", "ripenv", "fzmode", "kmpart", "nnhclu",
         "stacel", "mantst", "stmavg", "jtdecay", "jtsurf", "jtcal", "jtacc",
         "bjtc")
dsx <- c("fair", "tlvc", "tlvlr", "aorwt", "bdig", "bnft", "skp", "rptc",
         "rswp", "rshl", "zgll", "zgem", "mvtll", "mvtem", "tfrp", "lwss",
         "chmb", "gstr", "sptst")
spx <- c(spf, dsx)
