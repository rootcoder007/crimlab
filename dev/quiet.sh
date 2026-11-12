#!/bin/sh
# Runs a command only on a quiet machine: waits (checking every 60 s) until
# the 1-minute load average is at most 1.0, runs the command, then reads
# the load again. Exit status: the command's, or 3 when the load at the end
# was above 1.5 (the run overlapped other work: repeat it).
#   sh dev/quiet.sh CMD [ARG ...]
set -u
load() { cut -d ' ' -f 1 /proc/loadavg; }
until awk -v l="$(load)" 'BEGIN { exit !(l <= 1.0) }'; do
  echo "quiet: load $(load) at $(date +%H:%M), waiting"
  sleep 60
done
echo "quiet: start, load $(load) at $(date '+%Y-%m-%d %H:%M %Z')"
"$@"
rc=$?
end=$(load)
echo "quiet: end, load $end at $(date '+%Y-%m-%d %H:%M %Z')"
[ "$rc" -ne 0 ] && exit "$rc"
awk -v l="$end" 'BEGIN { exit !(l > 1.5) }' && exit 3
exit 0
