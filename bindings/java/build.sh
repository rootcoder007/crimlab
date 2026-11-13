#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Builds the Java binding with no build tool and no download: the shared
# library (cli/Makefile) and crimlab.jar (javac and jar, JDK 22 or later)
# under build/. With "test", runs the tests with the library given as
# -Dcrimlab.library; with "bundled", the jar also carries the library
# (native/<os>-<arch>/, as the Maven Central jar does) and the tests run
# with no library given, so they use that copy.
#   sh build.sh [test] [bundled]
set -eu
cd "$(dirname "$0")"
b=build
mkdir -p "$b"
nice -n 19 make -C ../../cli -j "${JOBS:-2}" all >"$b/make.log"
case $(uname -s) in
Darwin) lib=../../cli/libcrimlab.0.dylib os=macos ;;
*) lib=../../cli/libcrimlab.so.0 os=linux ;;
esac
f=$(java -XshowSettings:properties -version 2>&1 | grep -q 'os.name = Mac' && echo libcrimlab.dylib || echo libcrimlab.so)
rm -rf "$b/classes" "$b/tclasses" "$b/bclasses"
# shellcheck disable=SC2046 # no path has a space
javac --release 22 -Xlint:all -Werror -d "$b/classes" $(find src/main/java -name '*.java')
jar --create --file "$b/crimlab.jar" -C "$b/classes" .
# shellcheck disable=SC2046
javac --release 22 -Xlint:all -Werror -cp "$b/crimlab.jar" -d "$b/tclasses" $(find src/test/java -name '*.java')
for w in "$@"; do
  case $w in
  test)
    java --enable-native-access=ALL-UNNAMED -Dcrimlab.library="$lib" \
      -cp "$b/crimlab.jar:$b/tclasses" io.github.rootcoder007.crimlab.Tests ;;
  bundled)
    case $(uname -m) in
    x86_64 | amd64) arch=x86_64 ;;
    aarch64 | arm64) arch=arm64 ;;
    *) echo "build.sh: no bundled name for $(uname -m)" >&2; exit 2 ;;
    esac
    cp -R "$b/classes" "$b/bclasses"
    mkdir -p "$b/bclasses/native/$os-$arch"
    cp -L "$lib" "$b/bclasses/native/$os-$arch/$f"
    jar --create --file "$b/crimlab-bundled.jar" -C "$b/bclasses" .
    env -u CRIMLAB_LIBRARY java --enable-native-access=ALL-UNNAMED \
      -cp "$b/crimlab-bundled.jar:$b/tclasses" io.github.rootcoder007.crimlab.Tests ;;
  *) echo "build.sh: unknown step $w (test, bundled)" >&2; exit 2 ;;
  esac
done
echo "build.sh: $b/crimlab.jar"
