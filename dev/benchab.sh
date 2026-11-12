#!/bin/sh
# Old-versus-new timing, one function at a time, the two builds alternating
# so that load on the machine falls on both alike:
#   sh dev/benchab.sh OLDLIB NEWLIB
fns="vodtst orrisk outtst dspbm dyadrt gtmp cllr lrpost mixlr dbsrch cllrmin
lrcal tipplt eceplt kinidx patidx clrsmp expbnd misrd rtmhs incap incapr
sellab losbnd sntbnd hzdsel proof"
R_ENVIRON_USER=/dev/null
export R_ENVIRON_USER
for f in $fns; do
  o=$(nice Rscript dev/bench.R "$1" "$f")
  n=$(nice Rscript dev/bench.R "$2" "$f")
  echo "$o | $n"
done
