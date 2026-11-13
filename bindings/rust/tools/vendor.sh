#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Copies the kernels (every .cpp and .h of ../../src that does not include
# R's headers), cli/jrun.cpp with its headers, and the licence files into
# this crate, so that `cargo package` ships a crate that builds on its
# own. build.rs prefers ../../src when it is there; vendor/ is ignored by
# git. Rerun after any change to the kernels.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
top=$(cd "$here/../.." && pwd)
rm -rf "$here/vendor"
mkdir -p "$here/vendor/src" "$here/vendor/cli"
grep -L -E '#[[:space:]]*include[[:space:]]*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)' \
  "$top"/src/*.cpp "$top"/src/*.h | while read -r f; do cp "$f" "$here/vendor/src/"; done
cp "$top/cli/jrun.cpp" "$top/cli/jrun.h" "$top/cli/json.h" "$top/cli/thm.inc" "$here/vendor/cli/"
cp "$top/LICENSE" "$top/NOTICE" "$here/"
echo "vendor.sh: $(find "$here/vendor/src" -type f | wc -l | tr -d ' ') kernel files and cli/jrun.cpp"
