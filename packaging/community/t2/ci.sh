#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Builds the crimlab T2 package the way T2 builds packages, inside a T2
# userland (as root, in a container made from a T2 ISO's root file system):
# the T2 tree is the one mounted at /t2, the crimlab tree at /cl.
#   sh /cl/packaging/community/t2/ci.sh
# It checks the .desc with T2's Check-PkgFormat, fetches the release tarball
# from its [D] URL with T2's Download (which checks the [D] checksum), builds
# and installs the package with Emerge-Pkg, then again with SDECFG_DO_CHECK,
# so that the build runs the crimlab tests, and uses the installed command.
set -eu
s=/usr/src/t2-src
p=/cl/packaging/community/t2/package/scientific/crimlab

update-ca-certificates >/dev/null 2>&1 || true
cp -a /t2/. "$s/"
# the ExactBuild makefiles are an svn:externals of the T2 tree, which its
# git mirror lacks; T2's Config builds its Lua accelerator with them
[ -f "$s/misc/luabash/build/top.make" ] ||
  svn export -q https://svn.exactcode.de/exact-build/trunk/build@75 \
    "$s/misc/luabash/build"
mkdir -p "$s/package/scientific/crimlab"
cp "$p"/* "$s/package/scientific/crimlab/"
cd "$s"

echo "== Check-PkgFormat"
out=$(scripts/Check-PkgFormat crimlab)
if [ -n "$out" ]; then
  echo "$out"
  exit 1
fi

echo "== Config"
scripts/Config -oldconfig >/dev/null
scripts/Config PARALLEL="$(nproc)"
grep -q ' crimlab ' config/default/packages

echo "== Download (the [D] URL, not a mirror)"
scripts/Download -mirror none crimlab
ls -l download/mirror/c/crimlab-*

echo "== Emerge-Pkg crimlab"
gcc --version | head -n 1
scripts/Emerge-Pkg -optional-deps=no -parallel 1 crimlab
# the .cache this build wrote, against the shipped one
cat /var/adm/cache/crimlab
grep '^\[DEP\]' "$p/crimlab.cache" > /tmp/dep.old
grep '^\[DEP\]' /var/adm/cache/crimlab > /tmp/dep.new
diff /tmp/dep.old /tmp/dep.new || echo "the build's [DEP] lines differ from crimlab.cache"

echo "== Emerge-Pkg crimlab again, with its tests (SDECFG_DO_CHECK)"
scripts/Config DO_CHECK=1
scripts/Emerge-Pkg -optional-deps=no -parallel 1 -force crimlab
log=$(ls /var/adm/logs/*-crimlab.log)
grep -E 'checks, 0 failed' "$log"
head -n 3 /var/adm/packages/crimlab
grep -v ' var/adm/' /var/adm/flists/crimlab
# nothing but the manual page under man1 (the install wrapper once left an
# empty folder there)
test "$(find /usr/share/man/man1 -maxdepth 1 -name 'crimlab*' | wc -l)" -eq 1

echo "== use"
crimlab orrisk --or 2 --rate 0.1 --share 0.3 | grep 0.1475205
crimlab selftest -j | grep '"nfail":0'
echo "== t2 package ok"
