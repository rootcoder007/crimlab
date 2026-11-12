#!/bin/sh
# Runs a command only on a quiet machine: waits (checking every 60 s) until
# the 1-minute load average is at most QMAX (default 1.0), runs the
# command, then reads the load again. Exit status: the command's, or 3 when
# the load at the end was above QEND (default 1.5: the run overlapped other
# work, repeat it). On a machine whose idle load is higher (one that also
# runs interactive sessions) set both, and the scripts' own per-case rule
# (a case is timed again when the load rose by more than 1.5) does the work.
#   sh dev/quiet.sh CMD [ARG ...]
set -u
load() {
  if [ -r /proc/loadavg ]; then cut -d ' ' -f 1 /proc/loadavg
  else sysctl -n vm.loadavg | tr -d '{}' | awk '{ print $1 }'; fi
}
until awk -v l="$(load)" -v m="${QMAX:-1.0}" 'BEGIN { exit !(l <= m) }'; do
  echo "quiet: load $(load) at $(date +%H:%M), waiting"
  sleep 60
done
echo "quiet: start, load $(load) at $(date '+%Y-%m-%d %H:%M %Z')"
"$@"
rc=$?
end=$(load)
echo "quiet: end, load $end at $(date '+%Y-%m-%d %H:%M %Z')"
[ "$rc" -ne 0 ] && exit "$rc"
awk -v l="$end" -v m="${QEND:-1.5}" 'BEGIN { exit !(l > m) }' && exit 3
exit 0
