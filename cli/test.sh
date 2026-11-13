#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Tests of the command-line program:   sh test.sh ./crimlab
# test.tsv holds, per line, a function, its arguments as JSON and the
# result the Python package gives for them (numbers compared to a
# relative 1e-12, for other C libraries). Then the self-test of the
# kernels, the CSV and option input, and the errors.
# shellcheck disable=SC2016 # the backticks are part of the messages
set -u
c=${1:-./crimlab}
here=$(cd "$(dirname "$0")" && pwd)
n=0
bad=0
ck() { # name, got, want
  n=$((n + 1))
  if ! awk -v a="$2" -v b="$3" 'BEGIN {
      na = split(a, x, /[][{},:]/); nb = split(b, y, /[][{},:]/)
      if (na != nb) exit 1
      for (i = 1; i <= na; i++) {
        if (x[i] == y[i]) continue
        if (x[i] !~ /^-?[0-9.]+(e[-+]?[0-9]+)?$/ || y[i] !~ /^-?[0-9.]+(e[-+]?[0-9]+)?$/) exit 1
        d = x[i] - y[i]; if (d < 0) d = -d
        m = x[i] < 0 ? -x[i] : x[i]
        if (d > 1e-12 * m) exit 1
      }
    }'; then
    bad=$((bad + 1))
    printf 'FAIL %s\n  got:  %s\n  want: %s\n' "$1" "$2" "$3"
  fi
}
tab=$(printf '\t')
while IFS="$tab" read -r f j w; do
  ck "$f $j" "$(printf '%s' "$j" | "$c" "$f" -j -f -)" "$w"
done < "$here/test.tsv"
ck "selftest" "$("$c" selftest -j | sed 's/"n":[0-9]*,//')" '{"nfail":0,"report":""}'
ck "csv" "$(printf 'group,day,night\na,40,50\nb,60,50\n' | "$c" vodtst -j --ref b)" \
  '{"group":["a","b"],"ref":[false,true],"shday":[0.4,0.6],"shnight":[0.5,0.5],"ratio":[0.66666666666666663,1],"sign":[-1,1]}'
ck "csv short column" "$(printf 'done,open\n30,100\n45,\n' | "$c" losbnd -j)" '{"mean":37.5,"age":100,"lower":58.333333333333336,"bias":20.833333333333332}'
ck "options" "$("$c" losbnd -j --done 30,45,60,90,120 --open '[100, 150, 200]')" '{"mean":69,"age":150,"lower":99.375,"bias":30.375}'
ck "table" "$("$c" orrisk --or 2 --rate 0.1 --share 0.3 | sed -n 4p)" 'or                  2'
err() { # name, want message, exit status, command...
  w=$1 m=$2 s=$3
  shift 3
  got=$("$@" 2>&1 >/dev/null </dev/null)
  rc=$?
  ck "$w" "$rc $got" "$s $m"
}
err "unknown argument" 'orrisk: unknown argument `x`; the arguments are: or rate share.' 1 "$c" orrisk --or 2 --rate 0.1 --share 0.3 --x 1
err "missing argument" 'losbnd: `open` is missing.' 1 "$c" losbnd --done 1
err "check" 'orrisk: `rate` must be one number strictly between 0 and 1.' 1 "$c" orrisk --or 2 --rate 1 --share 0.3
err "no function" 'crimlab: no function nope (crimlab -l lists them).' 2 "$c" nope
err "bad JSON" 'crimlab: --p: the value is not JSON: [1,' 2 "$c" cllr --p '[1,'
"$c" -l | grep -q '^vodtst day night \[ref\] \[shift\]$' || { bad=$((bad + 1)); echo "FAIL -l"; }
"$c" -V | grep -q '^crimlab ' || { bad=$((bad + 1)); echo "FAIL -V"; }
# every function of the Python package has its entry in clab_jrun()
pyi="$here/../py/crimlab/__init__.py"
if [ -f "$pyi" ]; then
  n=$((n + 1))
  miss=$(sed -n '/^__all__ = \[/,/^]/p' "$pyi" | grep -o '"[a-z][a-z0-9]*"' | tr -d '"' |
    while read -r f; do "$c" -l | grep -q "^$f\( \|$\)" || printf '%s ' "$f"; done)
  [ -z "$miss" ] || { bad=$((bad + 1)); echo "FAIL not in clab_jrun(): $miss"; }
fi
echo "test.sh: $n checks, $bad failed"
[ "$bad" -eq 0 ]
