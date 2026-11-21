# The exported spatial (spf) and specialist (dsx) functions, whose shapes
# test-sp.R and test-spec.R check, the open data functions (odx), which
# test-odat.R checks on a copy of the service, the functions that ask the
# network (srx), and the placement volatility
# functions (cxx), which test-plcvol.R (plcvol, alrtcx), test-otpyr.R,
# test-cntmix.R, test-cxvol.R (cxvol, cxdml, cxyoy) and test-cxsim.R check;
# spx is all of them, left out of the shape test of test-crimlab.R.
spf <- c("hvdist", "eqproj", "eqinv", "cgstat", "cgmedn", "webpt", "sdell",
         "cvhull", "nnidx", "ripk", "ripenv", "fzmode", "kmpart", "nnhclu",
         "stacel", "mantst", "stmavg", "jtdecay", "jtsurf", "jtcal", "jtacc",
         "bjtc")
dsx <- c("fair", "tlvc", "tlvlr", "aorwt", "bdig", "bnft", "skp", "rptc",
         "rswp", "rshl", "zgll", "zgem", "mvtll", "mvtem", "tfrp", "lwss",
         "chmb", "gstr", "sptst")
spw <- c("spwt", "gearc", "gorg", "lmorn", "locgi", "dirmn")
sp3x <- c("hdbang", "crwalk", "rnnhcl")
odx <- c("odtab", "odls")
# the functions that ask a publisher's portal or the release page, which
# test-src.R and test-upd.R check on answers served to them
srx <- c("ckanres", "ckanls", "socres", "arcres", "srcget", "scres",
         "indtabs", "indiex", "ovrget", "updchk")
cxx <- c("plcvol", "alrtcx", "otpyr", "cntmix", "cxvol", "cxdml", "cxyoy",
         "cxsim")
spx <- c(spf, dsx, spw, sp3x, odx, srx, cxx)
