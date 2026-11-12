#!/bin/sh
# Mutation check of the self-test: each line of dev/mut.tsv is one
# deliberate bug in a scratch copy of src/core.cpp (kernel, text, its
# replacement; the text must occur exactly once). The standalone self-test
# is rebuilt against each and must fail. Run from the package root:
#   sh dev/mut.sh            prints one line per mutation and a summary
set -u
cxx=${CXX:-g++}
t=$(mktemp -d)
cp src/*.cpp src/*.h "$t/"
caught=0; missed=0
tab=$(printf '\t')
while IFS="$tab" read -r k from to; do
  case $k in ''|'#'*) continue ;; esac
  cp src/core.cpp "$t/core.cpp"
  F=$from T=$to perl -0pi -e '$n = s/\Q$ENV{F}\E/$ENV{T}/g;
    END { exit($n == 1 ? 0 : 3) }' "$t/core.cpp" ||
    { echo "MUT $k: text not found once"; missed=$((missed + 1)); continue; }
  if ! $cxx -std=c++17 -O1 -DCLAB_MAIN "$t/core.cpp" "$t"/t_*.cpp \
      "$t/selftest.cpp" -o "$t/st" 2>"$t/err"; then
    echo "MUT $k: does not compile"; missed=$((missed + 1)); continue
  fi
  if "$t/st" > "$t/out"; then
    echo "MUT $k: MISSED"; missed=$((missed + 1))
  else
    ids=$(grep ': want' "$t/out" | cut -d: -f1 | sort -u | head -5 | tr '\n' ' ')
    echo "MUT $k: caught, $(grep -o "[0-9]* checks, [0-9]* failed" "$t/out"); first ids: $ids"
    caught=$((caught + 1))
  fi
done < dev/mut.tsv
rm -rf "$t"
echo "MUT total: $caught caught, $missed missed"
[ "$missed" = 0 ]
