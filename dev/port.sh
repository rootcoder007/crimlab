#!/bin/sh
# The standalone self-test (no R) on many compilers, flags and machines.
# Run from the package root; the same commands run in .github/workflows/
# port.yml.
#   sh dev/port.sh flags          local compilers x flags (incl. sanitizers)
#   sh dev/port.sh mem            lowest ulimit -v, peak RSS, no allocation
#   sh dev/port.sh m32            native 32-bit x86, x87 and SSE2
#   sh dev/port.sh bench          ns per element of each kernel
#   sh dev/port.sh emu ARCH...    in a container under qemu: i386 i387
#                                 armv7 aarch64 riscv64 ppc64le s390x musl
#   sh dev/port.sh in             (inside a container) build and run
# ENGINE: podman (default) or docker. Prints one PORT line per run.
set -u
# every source that does not include R: the rule of cli/Makefile, so a new
# kernel is here the day it is written
src=$(grep -L -E '#[[:space:]]*include[[:space:]]*[<"](R\.h|Rinternals\.h|Rmath\.h|R_ext/|rchk\.h)' src/*.cpp | tr '\n' ' ')
# -Wno-psabi: on 32-bit ARM gcc notes, for every std::vector, that the
# calling convention changed in GCC 7.1; it is no warning about this code
base="-std=c++17 -fno-exceptions -fno-rtti -Wall -Wextra -pedantic -Wno-psabi -DCLAB_MAIN"
t=${TMPDIR:-/tmp}/clab-port.$$
mkdir -p "$t"
trap 'rm -rf "$t"' EXIT
rc=0

# one build and run: a label, the compiler, extra flags
one() {
  lab=$1; cc=$2; shift 2
  # shellcheck disable=SC2086
  if ! $cc $base "$@" $src -o "$t/st" -lm 2>"$t/w"; then
    echo "PORT $lab: build FAILED"; sed 's/^/  /' "$t/w" | head -20; rc=1
    return
  fi
  if [ -s "$t/w" ]; then
    echo "PORT $lab: warnings"; sed 's/^/  /' "$t/w" | head -20; rc=1
  fi
  if "$t/st" > "$t/o" 2>&1; then
    echo "PORT $lab: $(grep -o "[0-9]* checks, [0-9]* failed" "$t/o")"
  else
    echo "PORT $lab: FAILED $(grep -o "[0-9]* checks, [0-9]* failed" "$t/o")"; head -20 "$t/o"; rc=1
  fi
}

have() { command -v "$1" >/dev/null 2>&1; }

case ${1:-flags} in
flags)
  for cc in g++ clang++; do
    have $cc || continue
    for o in -O0 -O2 -Os; do one "$cc $o" $cc $o; done
    one "$cc -O2 fp-contract=off" $cc -O2 -ffp-contract=off
    one "$cc -O2 fp-contract=fast native" $cc -O2 -ffp-contract=fast \
      -march=native
    # only where the sanitizer runtime is installed
    if echo 'int main(){}' | $cc -x c++ -fsanitize=address,undefined - \
        -o "$t/san" 2>/dev/null; then
      # -Wno-maybe-uninitialized, for gcc (clang has no such option): with
      # the sanitizers on, gcc reports std::function inside its own <regex>
      # as possibly uninitialised
      mu=
      case $cc in g++*) mu=-Wno-maybe-uninitialized ;; esac
      # shellcheck disable=SC2086
      one "$cc asan+ubsan" $cc -O1 -g -fsanitize=address,undefined \
        -fno-sanitize-recover=all $mu
    else
      echo "PORT $cc asan+ubsan: skipped, no sanitizer runtime"
    fi
  done
  if gcc -std=c11 -Wall -Wextra -pedantic -fsyntax-only -x c src/crimlab.h
  then echo "PORT crimlab.h as C11: ok"; else echo "PORT crimlab.h: FAILED"; rc=1
  fi ;;
mem)
  cc=${CXX:-g++}
  # shellcheck disable=SC2086
  $cc $base -O2 $src -o "$t/st"
  echo "PORT peak RSS: $( (/usr/bin/time -f '%M KB' "$t/st" >/dev/null) 2>&1)"
  lo=0
  kb=1024
  while [ $kb -le 65536 ]; do
    if (ulimit -v $kb && "$t/st" >/dev/null 2>&1); then lo=$kb; break; fi
    kb=$((kb + 256))
  done
  echo "PORT lowest ulimit -v that passes: $lo KB"
  # what the kernels link against: libm and nothing that allocates
  $cc -std=c++17 -O2 -c src/core.cpp -o "$t/core.o"
  echo "PORT core.o undefined symbols: $(nm -u "$t/core.o" | awk '{print $2}' |
    tr '\n' ' ')"
  # shellcheck disable=SC2086
  $cc $base -DCLAB_WRAP -O2 $src -o "$t/stw" \
    -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc
  echo "PORT $("$t/stw" | grep allocations)" ;;
m32)
  for o in "-O2" "-O2 -mfpmath=387" "-O2 -msse2 -mfpmath=sse"; do
    # shellcheck disable=SC2086
    one "g++ -m32 $o" g++ -m32 $o
  done
  # the command and its whole test table with the arithmetic of x87, as a
  # 32-bit x86 distribution without SSE2 builds it (the package builders
  # for such systems ran it first, after a release was tagged)
  mkdir -p "$t/c32"
  cp -R cli src inst DESCRIPTION "$t/c32/"
  rm -rf "$t/c32/cli/build" "$t/c32/cli/crimlab" "$t/c32/cli"/libcrimlab.*
  if make -s -C "$t/c32/cli" -j4 check CXX="g++ -m32 -mfpmath=387" WERROR= > "$t/c32.log" 2>&1; then
    echo "PORT cli g++ -m32 x87: $(grep -o "[0-9]* checks, [0-9]* failed" "$t/c32.log")"
  else
    echo "PORT cli g++ -m32 x87: FAILED"; grep -A2 "^FAIL" "$t/c32.log" | cut -c1-300 | head -30
    tail -5 "$t/c32.log"; rc=1
  fi ;;
bench)
  cc=${CXX:-g++}
  # shellcheck disable=SC2086
  $cc $base -O2 $src -o "$t/st" && "$t/st" bench ;;
in)
  # inside a container: the compiler it has, plain and as small as can be
  # uname -m names the kernel's machine: a 32-bit x86 image on a 64-bit
  # host says x86_64, so the label adds the userland's word size
  cc=${CXX:-g++}
  m="$(uname -m)/$(getconf LONG_BIT)"
  one "$m $cc -O2" $cc -O2
  one "$m $cc -Os" $cc -Os
  [ -n "${XFLAGS:-}" ] && {
    # shellcheck disable=SC2086
    one "$m $cc -O2 $XFLAGS" $cc -O2 $XFLAGS
  } ;;
emu)
  shift
  eng=${ENGINE:-podman}
  for a in "$@"; do
    deb="apt-get update -qq && apt-get install -y -qq g++ >/dev/null"
    img=docker.io/library/debian:trixie-slim; x=""; pre=$deb
    case $a in
    i386) p=linux/386; x="-msse2 -mfpmath=sse" ;;
    i387) p=linux/386; x="-mfpmath=387" ;;
    armv7) p=linux/arm/v7 ;;
    # these have fused multiply-add; -std=c++17 leaves contraction off, so
    # one run turns it on
    aarch64) p=linux/arm64; x="-ffp-contract=fast" ;;
    riscv64) p=linux/riscv64; x="-ffp-contract=fast" ;;
    ppc64le) p=linux/ppc64le; x="-ffp-contract=fast" ;;
    s390x) p=linux/s390x; x="-ffp-contract=fast" ;;
    musl) p=linux/amd64; img=docker.io/library/alpine:3.22
      pre="apk add -q g++" ;;
    *) echo "PORT $a: unknown"; rc=1; continue ;;
    esac
    # label=disable: SELinux hosts refuse the bind mount otherwise
    $eng run --rm --platform $p --security-opt label=disable \
      -v "$PWD:/src:ro" -w /src -e XFLAGS="$x" \
      "$img" sh -c "$pre && sh dev/port.sh in" || rc=1
  done ;;
*) echo "usage: sh dev/port.sh flags|mem|bench|emu ARCH...|in"; rc=1 ;;
esac
exit $rc
