# v-bench: work in progress (stripped at merge)

Clone: l14 ~/work/crimlab-v, branch v-bench from main d9e160c.
Libraries: ~/tmp/benchlib (counterparts, bench, rmorie 1.4.3 from
r-universe, crimlab for the timed runs), ~/tmp/kt0 and ~/tmp/kt1
(instrumented builds of d9e160c and of the tree, dev/kcost.sh),
~/tmp/clab1 (plain build of the tree). Run R as R_ENVIRON_USER=/dev/null
R_PROFILE_USER=/dev/null. Gate: SRRLIB=~/tmp/Rlib.

## Done
- benchlib: bench, Familias, forrel, pedtools, comparison, sna, rmorie,
  relSim, DasGuptR, DemoDecomp, paramlink, tnet, netmem, NetworkToolbox,
  centiserve, idiolect (with its dependencies).
- tkrplot built against the Tcl/Tk 9.0.2 public headers from the upstream
  source tarballs (~/tmp/tkhdr/inc, via R_MAKEVARS_USER); no system package.
- forensim: not installable on l14. Its lazy-load attaches tcltk2 and
  tkrplot; without an X display tkrplot's Rplot_Init loops in Tcl 9's
  package search (stopped after 8 min; reproduced in 30 s with
  library(tcltk2); library(tkrplot)). No Xvfb on l14 (a system package).
- Counterparts checked by value: relSim lrSib/lrPC/lrMix, paramlink LR,
  idiolect:::Cllr/opt.llr, DasGuptR dgnpop, DemoDecomp horiuchi, sna, tnet,
  netmem (scripts ~/tmp/v1.R ... v8.R). dev/cpart.R and
  inst/bench/counterparts.csv written.
- Profile (callgrind, ~/tmp/prof): named 1e5 vodtst spends 83% of the
  .Call in name matching (text hashing); small calls spend ~4.7k
  instructions in seqs() (R's deferred conversion reads the print options).
- Binding: address-mode string sets, eager labels up to 64, incap labels
  written directly, labels symbol cached. Tests pass (1378); new
  tests/testthat/test-fast.R not yet run.
- dev/bench.R (run and merge modes, kinds, new counterparts), dev/bench-net.R
  (value checks, sna/tnet/netmem, merge), dev/benchrun.sh, dev/quiet.sh,
  README Speed prose, bench workflow packages.

## Next
1. cost before/after: sh dev/quiet.sh Rscript dev/cost.R ~/tmp/cost1.csv
   ~/tmp/kt1 ~/tmp/kt0 (log ~/tmp/cost1.log)
2. test bench.R and bench-net.R on one size, then install the tree into
   ~/tmp/benchlib and run: sh dev/benchrun.sh ~/tmp/benchlib 3 1800
   ~/tmp/benchruns
3. gate: SRRLIB=~/tmp/Rlib GATE="doc test lint cov check srr" sh dev/gate.sh
4. commits, bundle to /run/media/rootcoder/WD_BLACK/backups/crimlab/

## Open questions
- none
