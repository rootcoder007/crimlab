#!/bin/sh
# The speed tables: RUNS runs (default 3) of dev/bench.R and dev/bench-net.R
# on a quiet machine (dev/quiet.sh: started at a 1-minute load average of at
# most 1.0, repeated when it ends above 1.5), GAP seconds apart (default
# 1800, so that the runs are taken on different occasions), then merged
# into inst/bench/ and the README. Kept runs are in DIR (default
# bench-runs; a run already there is not repeated, so an interrupted
# series resumes).
#   sh dev/benchrun.sh LIB [RUNS] [GAP] [DIR]   (from the package root; LIB,
#                                                first on the library path,
#                                                holds crimlab and the
#                                                counterparts)
set -u
lib=$1
runs=${2:-3}
gap=${3:-1800}
dir=${4:-bench-runs}
R_ENVIRON_USER=/dev/null
R_PROFILE_USER=/dev/null
R_LIBS=$lib
# one process: no threads in BLAS or OpenMP code either
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
export R_ENVIRON_USER R_PROFILE_USER R_LIBS OMP_NUM_THREADS \
  OPENBLAS_NUM_THREADS MKL_NUM_THREADS
mkdir -p "$dir"
one() {
  # $1: script, $2: run directory, $3: the file a finished run leaves
  [ -f "$2/$3" ] && return 0
  until sh dev/quiet.sh Rscript "dev/$1" "$2" >"$2.log" 2>&1; do
    rc=$?
    if [ "$rc" -ne 3 ]; then
      echo "benchrun: $1 failed ($rc), see $2.log"
      exit "$rc"
    fi
    echo "benchrun: load above 1.5 at the end of $2, repeating"
    rm -rf "$2"
  done
  echo "benchrun: $2 kept"
}
k=1
while [ "$k" -le "$runs" ]; do
  one bench.R "$dir/bench$k" bench.csv
  one bench-net.R "$dir/net$k" net.csv
  if [ "$k" -lt "$runs" ] && [ ! -f "$dir/bench$((k + 1))/bench.csv" ]; then
    sleep "$gap"
  fi
  k=$((k + 1))
done
b=""
n=""
k=1
while [ "$k" -le "$runs" ]; do
  b="$b $dir/bench$k"
  n="$n $dir/net$k"
  k=$((k + 1))
done
# shellcheck disable=SC2086
Rscript dev/bench.R --merge $b && Rscript dev/bench-net.R --merge $n &&
  echo "benchrun: merged $runs runs"
