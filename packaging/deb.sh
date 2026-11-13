#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Builds the Debian packages (crimlab, libcrimlab0, libcrimlab-dev) from
# this checkout's HEAD: the source as git archive gives it, with
# packaging/debian as debian/, then dpkg-buildpackage. Needs git,
# build-essential and debhelper; the .deb files land in OUT (default: the
# current directory). With lintian installed, they are checked.
#   sh packaging/deb.sh [OUT]
set -eu
top=$(cd "$(dirname "$0")/.." && pwd)
out=$(mkdir -p "${1:-.}" && cd "${1:-.}" && pwd)
v=$(sed -n 's/^Version: *//p' "$top/DESCRIPTION")
w=$(mktemp -d)
trap 'rm -rf "$w"' EXIT
git -C "$top" archive --prefix="crimlab-$v/" -o "$w/crimlab_$v.orig.tar.gz" HEAD
tar -xzf "$w/crimlab_$v.orig.tar.gz" -C "$w"
cp -R "$top/packaging/debian" "$w/crimlab-$v/debian"
grep -q "^crimlab ($v-" "$w/crimlab-$v/debian/changelog" || {
  echo "deb.sh: debian/changelog has no entry for $v" >&2
  exit 1
}
(cd "$w/crimlab-$v" && dpkg-buildpackage -us -uc -b)
cp "$w"/*.deb "$out/"
if command -v lintian >/dev/null; then
  # the binary packages, not the .changes of an upload: these packages
  # go to no archive (so no ITP bug to close, and no distribution in the
  # changelog that Ubuntu's lintian would accept as its own)
  lintian --fail-on error,warning --suppress-tags initial-upload-closes-no-bugs "$w"/*.deb
fi
ls -l "$out"
