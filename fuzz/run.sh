#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Runs the fuzz target for a number of seconds from its seed corpus. A
# crash, a read or write out of bounds, undefined behaviour or a leak
# fails it, and the input is left in <dir>/run. An input that asks for
# hours of work (two thousand million iterations) or for more memory than
# there is does what it was asked: libFuzzer stops there (25 seconds,
# 2 GiB), the input is set aside in <dir>/run/slow and the run goes on.
#   fuzz/run.sh <dir> <seconds>
set -eu
d=${1:?usage: fuzz/run.sh <dir> <seconds>}
s=${2:?usage: fuzz/run.sh <dir> <seconds>}
mkdir -p "$d/run/slow" "$d/new"
cd "$d/run"
end=$(($(date +%s) + s))
while :; do
  left=$((end - $(date +%s)))
  [ "$left" -gt 0 ] || break
  rc=0
  # an allocation too large is std::bad_alloc, which clab_jrun reports
  ASAN_OPTIONS=allocator_may_return_null=1 ../fjrun -max_total_time="$left" -timeout=25 \
    -rss_limit_mb=2048 -max_len=4096 -print_final_stats=1 ../new ../corpus > log 2>&1 || rc=$?
  grep -a -E "^stat::|^SUMMARY|^Done" log || true
  [ "$rc" != 0 ] || break
  if ls crash-* leak-* > /dev/null 2>&1; then
    tail -n 60 log
    exit 1
  fi
  # neither a crash nor a slow input: the fuzzer itself failed
  n=0
  for f in timeout-* oom-*; do
    [ -f "$f" ] || continue
    mv "$f" slow/
    n=$((n + 1))
  done
  [ "$n" -gt 0 ] || { tail -n 60 log; exit 1; }
done
echo "fuzz/run.sh: no crash; $(find slow -type f | wc -l) slow or large inputs set aside"
