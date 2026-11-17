# The exported spatial (spf) and specialist (dsx) functions, whose shapes
# test-sp.R and test-spec.R check, and the open data function (odx), which
# test-odat.R checks on a copy of the service; spx is all of them, left out
# of the shape test of test-crimlab.R.
spf <- c("hvdist", "eqproj", "eqinv", "cgstat", "cgmedn", "webpt", "sdell",
         "cvhull", "nnidx", "ripk", "ripenv", "fzmode", "kmpart", "nnhclu",
         "stacel", "mantst", "stmavg", "jtdecay", "jtsurf", "jtcal", "jtacc",
         "bjtc")
dsx <- c("fair", "tlvc", "tlvlr", "aorwt", "bdig", "bnft", "skp", "rptc",
         "rswp", "rshl", "zgll", "zgem", "mvtll", "mvtem", "tfrp", "lwss",
         "chmb", "gstr", "sptst")
spw <- c("spwt", "gearc", "gorg", "lmorn", "locgi", "dirmn")
sp3x <- c("hdbang", "crwalk", "rnnhcl")
odx <- "odtab"
spx <- c(spf, dsx, spw, sp3x, odx)
