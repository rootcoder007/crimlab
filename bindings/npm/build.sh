#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# Builds dist/crimlab.wasm from the kernels (every .cpp of ../../src that
# does not include R's headers), cli/jrun.cpp and the command-line program
# cli/main.cpp (its main() exported as clab_cmain) with wasi-sdk
# (wasm32-wasip1), as a reactor: lib/wasm.js instantiates it and gives it
# the few WASI functions it calls. Run with WASI_SDK_PATH set to an
# unpacked wasi-sdk 34 (npm.yml downloads it and checks its SHA-256).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
want=34
sdk=${WASI_SDK_PATH:?build.sh: set WASI_SDK_PATH to wasi-sdk $want}
ver=$(sed -n 1p "$sdk/VERSION" 2>/dev/null | cut -d. -f1)
if [ "$ver" != "$want" ]; then
  echo "build.sh: wasi-sdk ${ver:-?} found, $want wanted" >&2
  exit 2
fi
v=$(sed -n 's/^Version: *//p' "$root/DESCRIPTION")
exp=$(printf ' -Wl,--export=%s' clab_jrun clab_jfns clab_cmain malloc free fflush)
ks=$(grep -L -E '#[[:space:]]*include[[:space:]]*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)' "$root"/src/*.cpp)
mkdir -p "$here/dist"
# shellcheck disable=SC2086 # the source list is split into words on purpose
"$sdk/bin/clang++" --target=wasm32-wasip1 --sysroot="$sdk/share/wasi-sysroot" \
  -mexec-model=reactor -std=c++17 -O2 -fno-exceptions \
  -Wall -Wextra -pedantic -Werror -ffile-prefix-map="$root"=. \
  -I"$root/src" -I"$root/cli" -DCLAB_VER="\"$v\"" -DCLAB_CMAIN \
  $ks "$root/cli/jrun.cpp" "$root/cli/main.cpp" \
  -Wl,-z,stack-size=1048576 -Wl,--strip-debug $exp \
  -o "$here/dist/crimlab.wasm"
cp "$root/LICENSE" "$root/NOTICE" "$here/"
ls -l "$here/dist"
