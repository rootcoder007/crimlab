#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# The release archive of the crimlab program for one platform, as
# release.yml builds it and install.sh installs it:
#   sh packaging/dist.sh TARGET OUT
# TARGET is linux-x86_64, linux-arm64, macos-x86_64 or macos-arm64, and
# must be this machine. Builds cli/ afresh (on Linux linked statically,
# so the program runs on any Linux of that CPU, glibc or musl), runs its
# tests, and writes OUT/crimlab-VERSION-TARGET.tar.gz, which holds
# crimlab-VERSION-TARGET/ with bin/crimlab, share/man/man1/crimlab.1,
# LICENSE and NOTICE. VER names another version (for dev/tabs.sh: the
# development version is not X.Y.Z, which install.sh asks for). Every file time is SOURCE_DATE_EPOCH (default: the
# commit time). GNU tar (TAR=gtar on macOS). Prints the archive's path.
# The Windows archive is built by release.yml with cl.
set -eu
t=${1:?usage: dist.sh TARGET OUT}
out=$(mkdir -p "${2:?usage: dist.sh TARGET OUT}" && cd "$2" && pwd)
top=$(cd "$(dirname "$0")/.." && pwd)
tar=${TAR:-tar}
v=${VER:-$(sed -n 's/^Version: *//p' "$top/DESCRIPTION")}
sde=${SOURCE_DATE_EPOCH:-$(git -C "$top" log -1 --format=%ct)}
case $(uname -s)-$(uname -m) in
  Linux-x86_64) me=linux-x86_64 ;;
  Linux-aarch64 | Linux-arm64) me=linux-arm64 ;;
  Darwin-x86_64) me=macos-x86_64 ;;
  Darwin-arm64) me=macos-arm64 ;;
  *) me=unknown ;;
esac
[ "$t" = "$me" ] || { echo "dist.sh: this machine is $me, not $t" >&2; exit 2; }
case $t in
  linux-*) lf=-static ;;
  *) lf= ; export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-11.0}" ;;
esac
make -C "$top/cli" clean >&2
make -C "$top/cli" -j"${JOBS:-4}" crimlab CXXFLAGS=-O2 LDFLAGS="$lf" >&2
strip "$top/cli/crimlab"
make -C "$top/cli" check >&2
n=crimlab-$v-$t
w=$(mktemp -d)
trap 'rm -rf "$w"' EXIT
mkdir -p "$w/$n/bin" "$w/$n/share/man/man1"
cp "$top/cli/crimlab" "$w/$n/bin/"
cp "$top/cli/crimlab.1" "$w/$n/share/man/man1/"
cp "$top/LICENSE" "$top/NOTICE" "$w/$n/"
(cd "$w" && "$tar" --sort=name --owner=0 --group=0 --numeric-owner \
  --mtime="@$sde" -cf - "$n") | gzip -n -9 > "$out/$n.tar.gz"
echo "$out/$n.tar.gz"
