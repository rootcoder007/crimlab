#!/bin/sh
# Old-versus-new timing, one function at a time, the two builds alternating
# so that load on the machine falls on both alike:
#   sh dev/benchab.sh OLDLIB NEWLIB
fns=$(cut -f1 dev/bench.txt)
R_ENVIRON_USER=/dev/null
export R_ENVIRON_USER
for f in $fns; do
  o=$(nice Rscript dev/bench.R "$1" "$f")
  n=$(nice Rscript dev/bench.R "$2" "$f")
  echo "$o | $n"
done
