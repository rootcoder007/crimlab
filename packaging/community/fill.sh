#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Writes the community recipes of one release. In this tree the recipes
# carry the version 0.0.0 and checksums of zeros, so that no commit has to
# follow a release; this script fills them in from the archives of the
# release, and the result is what is sent to each outside repository.
#   sh packaging/community/fill.sh VERSION DIR OUT
# DIR holds crimlab-VERSION.tar.gz and crimlab-VERSION-windows-x86_64.zip;
# OUT is made, a copy of this directory with every recipe filled in.
# Checksums: SHA-256 and SHA-512 of the source tarball, T2's SHA-224 of the
# decompressed tarball, and for winget the SHA-256 of the Windows archive.
set -eu
[ $# -eq 3 ] || { echo "usage: fill.sh VERSION DIR OUT" >&2; exit 2; }
v=$1
d=$2
o=$3
here=$(cd "$(dirname "$0")" && pwd)
src=$d/crimlab-$v.tar.gz
win=$d/crimlab-$v-windows-x86_64.zip
for f in "$src" "$win"; do
  [ -f "$f" ] || { echo "fill.sh: $f is not there" >&2; exit 2; }
done
case $v in
  [0-9]*.[0-9]*.[0-9]*) ;;
  *) echo "fill.sh: $v is not a version" >&2; exit 2 ;;
esac

s256=$(sha256sum < "$src" | cut -d' ' -f1)
s512=$(sha512sum < "$src" | cut -d' ' -f1)
s224=$(gzip -dc "$src" | sha224sum | cut -d' ' -f1)
w256=$(sha256sum < "$win" | cut -d' ' -f1 | tr 'a-f' 'A-F')
z56=$(printf '%056d' 0)
z64=$(printf '%064d' 0)
z128=$(printf '%0128d' 0)

mkdir -p "$o"
cp -R "$here"/. "$o"/
rm -f "$o/fill.sh"
wg=$o/winget/manifests/r/rootcoder007/crimlab
mv "$wg/0.0.0" "$wg/$v"

# the longest run of zeros first: each shorter one is part of it. The
# version is 0.0.0 with no digit or dot before it and no digit after it,
# nor a dot and a digit (crimlab-0.0.0.tar.gz is one, 0.0.0.0 is not);
# twice, for two versions that share a character between them
ver="s/(^|[^0-9.])0\\.0\\.0(\\.?([^0-9.]|\$))/\\1$v\\2/g"
find "$o" -type f | while read -r f; do
  case $f in
    "$wg"/*) h=$w256 ;;
    *) h=$s256 ;;
  esac
  sed -E -e "s/$z128/$s512/g" -e "s/$z64/$h/g" -e "s/$z56/$s224/g" \
    -e "$ver" -e "$ver" "$f" > "$f.new"
  # written back into the file, which keeps its mode (t2/ci.sh runs)
  cat "$f.new" > "$f"
  rm -f "$f.new"
done

# nothing is left unfilled
if grep -rnE "(^|[^0-9.])0\\.0\\.0(\\.?([^0-9.]|\$))|$z56" "$o"; then
  echo "fill.sh: the lines above were not filled in" >&2
  exit 1
fi
echo "fill.sh: recipes of crimlab $v in $o"
