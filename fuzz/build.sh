#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Builds the fuzz target fuzz/fjrun.cpp with libFuzzer, AddressSanitizer and
# UndefinedBehaviorSanitizer (clang), and its seed corpus from the
# command-line test table: one file per row, the function name on the
# first line and its JSON arguments after it.
#   fuzz/build.sh <dir>
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
out=${1:?usage: fuzz/build.sh <dir>}
mkdir -p "$out/corpus"
out=$(cd "$out" && pwd)
san="-fsanitize=fuzzer-no-link,address,undefined -fno-sanitize-recover=undefined"
make -C "$here/cli" -j"${CLJ:-4}" B="$out/o" libcrimlab.a CXX="${CXX:-clang++}" \
  CXXFLAGS="-O1 -g -fno-omit-frame-pointer $san" WERROR= > "$out/build.log" 2>&1 ||
  { tail -20 "$out/build.log" >&2; exit 1; }
mv "$here/cli/libcrimlab.a" "$out/libcrimlab.a"
"${CXX:-clang++}" -std=c++17 -O1 -g -fno-omit-frame-pointer \
  -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=undefined \
  -I"$here/cli" -I"$here/src" "$here/fuzz/fjrun.cpp" "$out/libcrimlab.a" -o "$out/fjrun"
awk -F'\t' -v d="$out/corpus" '
  NF >= 2 { f = sprintf("%s/%05d", d, NR); printf "%s\n%s", $1, $2 > f; close(f) }
' "$here/cli/test.tsv"
echo "fuzz/build.sh: $out/fjrun, $(ls "$out/corpus" | wc -l) seeds"
