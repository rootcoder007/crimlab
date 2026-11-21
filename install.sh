#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# Installs the crimlab command and its manual page from a GitHub release
# of rootcoder007/crimlab, on Linux and macOS, into ~/.local (no sudo).
#
#   curl -fsSL https://github.com/rootcoder007/crimlab/releases/latest/download/install.sh | sh
#   sh install.sh [--prefix DIR] [--version X.Y.Z]
#
# In order, stopping at the first failure:
#   1. picks the archive for this OS and CPU (linux or macos, x86_64 or
#      arm64);
#   2. downloads the release's SHA256SUMS (and, with CRIMLAB_SIGNER set,
#      checks its SSH signature SHA256SUMS.sig);
#   3. downloads the archive and checks its SHA-256 against SHA256SUMS;
#   4. checks every path in the archive stays inside its top directory;
#   5. copies bin/ and share/ into the prefix.
# Nothing from the archive runs before step 3 passes. POSIX sh: dash,
# busybox and bash all run it.
#
# Environment (the options win):
#   CRIMLAB_PREFIX   where to install (default: $HOME/.local)
#   CRIMLAB_VERSION  a release such as 0.1.0 (default: the latest)
#   CRIMLAB_SIGNER   another release key, as an allowed-signers line
#                    ("crimlab ssh-ed25519 AAAA..."), in place of the
#                    one below (clabsigner, from SECURITY.md), which
#                    SHA256SUMS.sig must verify against
#   CRIMLAB_VERIFY   0: go on without checking who signed the release
#   CRIMLAB_BASE     the releases URL (default:
#                    https://github.com/rootcoder007/crimlab/releases);
#                    https, or http to localhost for tests
#
# The whole script is one function called on its last line, so a
# download cut short runs nothing.
set -eu

clabrepo=rootcoder007/crimlab
# the release signing key (SECURITY.md)
clabsigner='crimlab ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDCgdm/W6LtsPET14oTqb1lOFga52VG4HSJTtQHfLTFc'

# On a terminal the steps are numbered, results are marked and coloured,
# the archive shows a progress bar and a table sums the files up; piped
# or in CI (no terminal, TERM=dumb or NO_COLOR set) every line is plain
# "crimlab install: ..." text, as it always was.
clabtty=0
if [ -t 2 ] && [ "${TERM:-dumb}" != dumb ] && [ -z "${NO_COLOR:-}" ]; then
  clabtty=1
fi
clabstep=0
say() {
  if [ "$clabtty" = 0 ]; then
    printf 'crimlab install: %s\n' "$*" >&2
    return
  fi
  case $* in
    error:*) printf '\033[1;31m  error\033[0m %s\n' "${*#error: }" >&2 ;;
    *"not checked"*) printf '\033[33m   note\033[0m %s\n' "$*" >&2 ;;
    *": good"* | *" matches "*) printf '\033[32m     ok\033[0m %s\n' "$*" >&2 ;;
    *) printf '        %s\n' "$*" >&2 ;;
  esac
}
# step <title> [plain line]: one of the five steps; the plain line is
# what is said where there is no terminal
step() {
  clabstep=$((clabstep + 1))
  if [ "$clabtty" = 1 ]; then
    printf '\033[1m[%d/5]\033[0m %s\n' "$clabstep" "$1" >&2
  else
    [ -z "${2:-}" ] || say "$2"
  fi
}
# the size of a file as people read it
hsize() {
  wc -c < "$1" | awk '{ b = $1
    if (b >= 1048576) printf "%.1f MiB", b / 1048576
    else if (b >= 1024) printf "%.1f KiB", b / 1024
    else printf "%d B", b }'
}
die() { say "error: $*"; say "nothing was installed"; exit 1; }

fetch() { # fetch URL FILE [bar]: with bar, a progress bar on a terminal
  clabq=-sS
  clabw=-q
  if [ "${3:-}" = bar ] && [ "$clabtty" = 1 ]; then
    clabq=-#
    clabw='-q --show-progress'
  fi
  if command -v curl >/dev/null 2>&1; then
    case $1 in
      https://*) curl -fL "$clabq" --proto '=https' --tlsv1.2 -o "$2" "$1" ;;
      *) curl -f "$clabq" --proto '=http' -o "$2" "$1" ;;
    esac
  elif command -v wget >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    wget $clabw -O "$2" "$1"
  else
    die "neither curl nor wget is installed"
  fi
}

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v openssl >/dev/null 2>&1; then
    openssl dgst -sha256 -r "$1" | awk '{print $1}'
  else
    die "no sha256sum, shasum or openssl to check the download with"
  fi
}

clabmain() {
  prefix=${CRIMLAB_PREFIX:-}
  ver=${CRIMLAB_VERSION:-}
  signer=${CRIMLAB_SIGNER:-$clabsigner}
  # CRIMLAB_VERIFY=0: go on without checking who signed the release
  [ "${CRIMLAB_VERIFY:-1}" != 0 ] || signer=
  base=${CRIMLAB_BASE:-https://github.com/$clabrepo/releases}
  while [ $# -gt 0 ]; do
    case $1 in
      --prefix) [ $# -ge 2 ] || die "--prefix needs a directory"; prefix=$2; shift 2 ;;
      --prefix=*) prefix=${1#--prefix=}; shift ;;
      --version) [ $# -ge 2 ] || die "--version needs a version"; ver=$2; shift 2 ;;
      --version=*) ver=${1#--version=}; shift ;;
      -h | --help) say "usage: sh install.sh [--prefix DIR] [--version X.Y.Z]"; exit 0 ;;
      *) die "unknown option: $1" ;;
    esac
  done
  if [ -z "$prefix" ]; then
    [ -n "${HOME:-}" ] || die "HOME is not set; give --prefix"
    prefix=$HOME/.local
  fi
  ver=${ver#v}
  vre='[0-9]+\.[0-9]+\.[0-9]+'
  if [ -n "$ver" ] && ! printf '%s\n' "$ver" | grep -Eqx "$vre"; then
    die "not a release version: $ver"
  fi
  base=${base%/}
  case $base in
    https://*) ;;
    *) printf '%s\n' "$base" |
         grep -Eq '^http://(127\.0\.0\.1|localhost)(:[0-9]{1,5})?(/[^@]*)?$' ||
         die "CRIMLAB_BASE must be https (or http to localhost): $base" ;;
  esac

  # 1. the archive for this machine
  case $(uname -s) in
    Linux) os=linux ;;
    Darwin) os=macos ;;
    MINGW* | MSYS* | CYGWIN*) die "on Windows use install.ps1" ;;
    *) die "no release archive for $(uname -s); build from source (cli/)" ;;
  esac
  case $(uname -m) in
    x86_64 | amd64) arch=x86_64 ;;
    aarch64 | arm64) arch=arm64 ;;
    *) die "no release archive for the $(uname -m) CPU; build from source (cli/)" ;;
  esac
  # an x86-64 shell under Rosetta on Apple silicon: take the native build
  if [ "$os" = macos ] && [ "$arch" = x86_64 ] &&
     [ "$(sysctl -n sysctl.proc_translated 2>/dev/null || echo 0)" = 1 ]; then
    arch=arm64
  fi
  target=$os-$arch
  if [ -n "$ver" ]; then url=$base/download/v$ver; else url=$base/latest/download; fi

  tmp=$(mktemp -d 2>/dev/null || mktemp -d -t crimlab)
  if [ -z "$tmp" ] || [ ! -d "$tmp" ]; then die "cannot make a temporary directory"; fi
  trap 'rm -rf "$tmp"' EXIT
  trap 'exit 1' INT TERM HUP

  # 2. the sums
  step "Fetch SHA256SUMS" "fetching $url/SHA256SUMS"
  fetch "$url/SHA256SUMS" "$tmp/SHA256SUMS" || die "cannot download SHA256SUMS"
  names=$(awk '{f=$2; sub(/^\*/, "", f); print f}' "$tmp/SHA256SUMS" |
    grep -Ex "crimlab-$vre-$target\.tar\.gz" || true)
  if [ -n "$ver" ]; then
    name=crimlab-$ver-$target.tar.gz
    printf '%s\n' "$names" | grep -qxF "$name" || die "SHA256SUMS of v$ver lists no $name"
  else
    [ -n "$names" ] || die "the latest release lists no archive for $target"
    [ "$(printf '%s\n' "$names" | wc -l | tr -d ' ')" = 1 ] ||
      die "SHA256SUMS lists more than one archive for $target"
    name=$names
    ver=${name#crimlab-}
    ver=${ver%-"$target".tar.gz}
  fi
  want=$(awk -v n="$name" '{f=$2; sub(/^\*/, "", f); if (f == n) print tolower($1)}' \
    "$tmp/SHA256SUMS")
  printf '%s\n' "$want" | grep -Eqx '[0-9a-f]{64}' ||
    die "SHA256SUMS has no single valid SHA-256 for $name"
  step "Check who signed it"
  if [ -n "$signer" ]; then
    command -v ssh-keygen >/dev/null 2>&1 ||
      die "ssh-keygen (OpenSSH) is not installed, so who signed this release cannot be checked; install it, or set CRIMLAB_VERIFY=0 to go on without"
    fetch "$url/SHA256SUMS.sig" "$tmp/SHA256SUMS.sig" || die "cannot download SHA256SUMS.sig"
    printf '%s\n' "$signer" > "$tmp/signers"
    ssh-keygen -Y verify -f "$tmp/signers" -I crimlab -n file \
      -s "$tmp/SHA256SUMS.sig" < "$tmp/SHA256SUMS" > "$tmp/ssh.log" 2>&1 ||
      { cat "$tmp/ssh.log" >&2; die "the SSH signature of SHA256SUMS does not verify"; }
    say "SSH signature of SHA256SUMS: good"
  else
    say "CRIMLAB_VERIFY=0: who signed this release is not checked"
  fi

  # 3. the archive, checked before anything else touches it
  step "Fetch $name" "fetching $url/$name"
  fetch "$url/$name" "$tmp/$name" bar || die "cannot download $name"
  got=$(sha256 "$tmp/$name")
  [ "$got" = "$want" ] || die "SHA-256 mismatch for $name: SHA256SUMS says $want, the download is $got"
  say "SHA-256 of $name matches SHA256SUMS"

  step "Check and unpack"
  # 4. every entry inside crimlab-<version>-<target>/
  top=crimlab-$ver-$target
  tar -tzf "$tmp/$name" > "$tmp/list" || die "cannot read $name"
  ! grep -Evq "^$top(/|\$)" "$tmp/list" || die "$name has entries outside $top/"
  ! grep -Eq '(^|/)\.\.(/|$)' "$tmp/list" || die "$name has a .. path"
  mkdir "$tmp/x"
  tar -xzf "$tmp/$name" -C "$tmp/x" || die "cannot unpack $name"
  [ -f "$tmp/x/$top/bin/crimlab" ] || die "$name has no bin/crimlab"

  step "Install into $prefix"
  # 5. into the prefix, the command last
  { mkdir -p "$prefix" 2>/dev/null && [ -w "$prefix" ]; } ||
    die "cannot write to $prefix; choose another with --prefix (this script never uses sudo)"
  for d in share bin; do
    mkdir -p "$prefix/$d"
    cp -R "$tmp/x/$top/$d/." "$prefix/$d/" || die "copy into $prefix/$d failed (it may be incomplete)"
  done
  if [ "$clabtty" = 1 ]; then
    printf '\n\033[1m %-40s %-10s %s\033[0m\n' File Size Check >&2
    printf ' %-40s %-10s %s\n' SHA256SUMS "$(hsize "$tmp/SHA256SUMS")" \
      "signature: see step 2" "$name" "$(hsize "$tmp/$name")" \
      "SHA-256 matches" >&2
    printf '\n\033[1;32m crimlab %s\033[0m (%s) installed into %s\n' \
      "$ver" "$target" "$prefix" >&2
  else
    say "installed crimlab $ver ($target) into $prefix"
  fi
  "$prefix/bin/crimlab" -V >&2 || die "the installed crimlab does not run"
  case ":${PATH:-}:" in
    *":$prefix/bin:"*) ;;
    *) say "add $prefix/bin to PATH:  export PATH=\"$prefix/bin:\$PATH\"" ;;
  esac
}

clabmain "$@"
