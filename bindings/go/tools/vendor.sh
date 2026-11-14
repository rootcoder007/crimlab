#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Copies the kernels (every .cpp and .h of ../../src that does not include
# R's headers, the rule py/setup.py follows) and cli/jrun.cpp with its
# headers into csrc/, writes one k_<file>.cpp per source that includes
# it, and writes fns.go, one Go function per entry of the table in
# jrun.cpp. cgo compiles only the files in the package's own folder, and
# a module fetched with `go get` holds only bindings/go, so the copy is
# committed; pkg.yml fails when it is out of date. Rerun after any change
# to the kernels or to cli/jrun.cpp.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
top=$(cd "$here/../.." && pwd)
rm -rf "$here/csrc"
rm -f "$here"/k_*.cpp
mkdir -p "$here/csrc/src" "$here/csrc/cli"
grep -L -E '#[[:space:]]*include[[:space:]]*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)' \
  "$top"/src/*.cpp "$top"/src/*.h | while read -r f; do cp "$f" "$here/csrc/src/"; done
cp "$top/cli/jrun.cpp" "$top/cli/jrun.h" "$top/cli/json.h" "$top"/cli/*.inc "$here/csrc/cli/"
for f in "$here"/csrc/src/*.cpp "$here/csrc/cli/jrun.cpp"; do
  d=$(basename "$(dirname "$f")")
  b=$(basename "$f" .cpp)
  printf '// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception\n// Made by tools/vendor.sh; do not edit.\n#include "csrc/%s/%s.cpp"\n' "$d" "$b" > "$here/k_$b.cpp"
done
{
  printf '// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception\n\n'
  printf '// Made by tools/vendor.sh from the table of cli/jrun.cpp; do not edit.\n\n'
  printf 'package crimlab\n'
  awk -F'"' '
    /const Ent ents\[\] = \{/ { on = 1; next }
    on && /^};/ { exit }
    on && NF >= 5 {
      n = $2; a = $4
      d = a == "" ? "no arguments" : "arguments: " a
      printf "\n// %s runs %s() of the R package (%s; ? marks an optional one).\n", toupper(substr(n, 1, 1)) substr(n, 2), n, d
      printf "func %s(args string) (string, error) { return Jrun(\"%s\", args) }\n", toupper(substr(n, 1, 1)) substr(n, 2), n
    }' "$top/cli/jrun.cpp"
} > "$here/fns.go"
cp "$top/LICENSE" "$top/NOTICE" "$here/"
echo "vendor.sh: $(find "$here/csrc/src" -type f | wc -l | tr -d ' ') kernel files, $(grep -c '^func' "$here/fns.go") functions"
