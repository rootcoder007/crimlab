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

## Change of machine (coordinator, 17:56)
Timings of record are taken on the Mac (Apple M2, 8 GB, Homebrew R 4.6.1),
not on l14; builds, tests and the gate stay on l14. Mac tree:
/Volumes/VSR/rootcoderfiles/crimlab-bench/v/crimlab (clone of the l14
bundle), libraries v/lib (crimlab and my counterparts), v/kt0, v/kt1, the
shared read-only lib crimlab-bench/lib. One timing job at a time across
agents: v/locked.sh takes crimlab-bench/lock (mkdir), writes owner, removes
it on exit. The Mac idles at a load of about 2 (the sessions run there):
each case records its load before and after and is timed again when it
rose by more than 1.5. No l14 timing run was completed (the cost runs were
stopped mid-way): no l14 numbers are kept.

## State at hand-over (2026-11-12 23:15)
- Tables of record taken on the Mac 21:38, 22:18, 22:58 (each under the
  lock, about 10 minutes); merged into inst/bench/bench.csv, bench.md,
  net.csv and the README lists. Rejected and kept out: an l14 cost run
  (stopped), the first Mac bench/net runs (one stopped by me at 19:22 to
  let the spatial table have the lock, one with forensim loaded, which
  opens a Tk window on the desktop; one net run that swapped in sna's
  betweenness at 1e4 and one whose eigenvector and modularity checks
  were too strict).
- Gate on l14 (rebased on abbe4f4): test 2921/0/0/0, lint 0, cov 100,
  gp 0 failed, srr TODO 0, check 2 NOTEs; dev/port.sh flags all 0 failed
  (g++ asan skipped: no sanitizer runtime); R 4.1.3 container builds and
  loads.

## Next
- nothing in this unit; merge strips this file.

## Open questions
- none
