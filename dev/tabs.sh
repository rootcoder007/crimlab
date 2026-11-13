#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Runs the install command of every channel the website shows, each in a
# clean container, against artifacts built here from this working tree
# (until the first release the registries do not have crimlab, so each
# command points at a local stand-in: the release archives served on
# localhost for install.sh, a folder of wheels for pip, the packed
# tarball for npm, the packaged crate for cargo, a file Maven repository,
# a file Go module proxy, a local tap for Homebrew, a file CRAN-like
# repository for R, the built .deb and .rpm files for apt and dnf). Each
# runs `orrisk --or 2 --rate 0.1 --share 0.3` (or its call in that
# language) and must print exposed = 0.1475205. One PASS, FAIL or SKIP
# line per channel; exit status 1 when any failed.
#   sh dev/tabs.sh [CHANNEL ...]
# Channels: git sh busybox ps1 docker brew scoop pip npm cargo go maven r
# apt dnf (default: all). Needs podman (or CTR=docker), python3 (the
# local server and the Go module proxy), tar and sha256sum. Works in
# ${TMPDIR:-/tmp}/clabtabs (logs: DIR/log/CHANNEL.log). The release
# archive is built for this machine's CPU.
set -u
top=$(cd "$(dirname "$0")/.." && pwd)
D=${TMPDIR:-/tmp}/clabtabs
ctr=${CTR:-podman}
command -v "$ctr" >/dev/null 2>&1 || ctr=docker
V=0.0.0 # the release archives' version: install.sh asks for X.Y.Z
P=${TABSPORT:-8471}
B=http://127.0.0.1:$P
K='0\.14752(05|049)'
WASI=https://github.com/WebAssembly/wasi-sdk/releases/download/wasi-sdk-34/wasi-sdk-34.0-arm64-linux.tar.gz
WASISUM=f7e243dff54d60bcc576e94d6166b69f410f2500ae4a9ceef34315be10e77971
chans=${*:-git sh busybox ps1 docker brew scoop pip npm cargo go maven r apt dnf}
PWSH=https://github.com/PowerShell/PowerShell/releases/download/v7.7.0/powershell-7.7.0-linux-arm64.tar.gz
PWSHSUM=e59c89776eaa59790d44b3edd14277435225ab08b6e7603763101b89de83d935
case $(uname -m) in
  x86_64 | amd64)
    arch=x86_64 plat=amd64
    WASI=$(printf '%s' "$WASI" | sed 's/arm64-linux/x86_64-linux/')
    WASISUM=b761e3a0721dbae9c09a0059e5fdb2bf917d1b4a8a7b430fb3b5aafb0984b2c4
    PWSH=$(printf '%s' "$PWSH" | sed 's/linux-arm64/linux-x64/')
    PWSHSUM=c89d9c682d2c403dfe726f3df948e359675235a5bd6fab535776de69729ca0cb ;;
  *) arch=arm64 plat=arm64 ;;
esac
# get URL SHA256 FILE: a download into the cache, kept while it matches
get() {
  if ! printf '%s  %s\n' "$2" "$D/cache/$3" | sha256sum -c - >/dev/null 2>&1; then
    curl -fsSL -o "$D/cache/$3" "$1" || return 1
    printf '%s  %s\n' "$2" "$D/cache/$3" | sha256sum -c - >/dev/null 2>&1 || return 1
  fi
}

# what a container wrote as root, removed by a container
rm -rf "$D/src" "$D/art" "$D/rel" "$D/log" 2>/dev/null ||
  "$ctr" run --rm -v "$D:/t" docker.io/library/alpine:3.22 rm -rf /t/src /t/art /t/rel /t/log
mkdir -p "$D/src" "$D/art" "$D/rel/download/v$V" "$D/log" "$D/cache"
bad=0

# a container: run NAME IMAGE SCRIPT (the script on standard input),
# /t is $D; the host network, for the server on localhost; this
# machine's CPU (an image pulled for another one would run emulated)
run() {
  timeout 1800 "$ctr" run --rm -i --network host --platform "linux/$plat" -v "$D:/t" \
    "$2" sh > "$D/log/$1.log" 2>&1 <<X
set -e
$3
X
}
res() { # res CHANNEL STATUS: PASS when the log shows the result
  if [ "$2" = 0 ] && grep -Eq "$K" "$D/log/$1.log"; then
    echo "PASS $1"
  else
    bad=1
    echo "FAIL $1 (rc=$2): $(tail -3 "$D/log/$1.log" | tr '\n' ' ' | cut -c1-240)"
  fi
}
want() { case " $chans " in *" $1 "*) return 0 ;; esac; return 1; }

# the working tree as a fresh git repository (deb.sh and the COPR
# Makefile take git archive of HEAD)
(cd "$top" && git ls-files -co --exclude-standard -z | xargs -0 tar -cf - 2>/dev/null) |
  tar -xf - -C "$D/src"
git -C "$D/src" init -q && git -C "$D/src" add -A &&
  git -C "$D/src" -c user.name=tabs -c user.email=tabs@localhost commit -qm tabs
sde=$(git -C "$top" log -1 --format=%ct)

# the release: the archive for this CPU, the source tarball, the scripts,
# SHA256SUMS; served on localhost
if want sh || want busybox || want brew; then
  run build-rel docker.io/library/debian:trixie "apt-get update -qq; apt-get install -y -qq g++ make git >/dev/null
cp -R /t/src /tmp/s; cd /tmp/s
VER=$V SOURCE_DATE_EPOCH=$sde sh packaging/dist.sh linux-$arch /t/rel/download/v$V" ||
    echo "build-rel failed: $(tail -3 "$D/log/build-rel.log" | tr '\n' ' ')"
  (cd "$D/src" && git archive --prefix="crimlab-$V/" HEAD | gzip -n -9 > "$D/rel/download/v$V/crimlab-$V.tar.gz")
  cp "$top/install.sh" "$top/install.ps1" "$D/rel/download/v$V/"
  (cd "$D/rel/download/v$V" && sha256sum -- * > SHA256SUMS)
  mkdir -p "$D/rel/latest" && cp -R "$D/rel/download/v$V" "$D/rel/latest/download"
fi
(cd "$D/rel" && exec python3 -m http.server "$P" --bind 127.0.0.1) > "$D/log/server.log" 2>&1 &
srv=$!
trap 'kill "$srv" 2>/dev/null' EXIT
trap 'exit 1' INT TERM HUP
sleep 1

if want git; then
  run git docker.io/library/debian:trixie-slim "apt-get update -qq; apt-get install -y -qq git g++ make ca-certificates >/dev/null
git config --global --add safe.directory '*'
git clone -q /t/src crimlab
make -C crimlab/cli -j4 crimlab >/dev/null
crimlab/cli/crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res git $?
fi

if want sh; then
  # debian's /bin/sh is dash
  run sh docker.io/library/debian:trixie-slim "apt-get update -qq; apt-get install -y -qq curl ca-certificates >/dev/null
ls -l /bin/sh
curl -fsSL $B/latest/download/install.sh | CRIMLAB_BASE=$B sh
~/.local/bin/crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res sh $?
fi

if want busybox; then
  run busybox docker.io/library/alpine:3.22 "wget -qO- $B/latest/download/install.sh | CRIMLAB_BASE=$B sh
~/.local/bin/crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res busybox $?
fi

if want ps1; then
  # Windows only: here PowerShell parses it, nothing more
  if get "$PWSH" "$PWSHSUM" pwsh.tgz; then
    run ps1 docker.io/library/debian:trixie-slim "export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
mkdir /tmp/p; tar -xzf /t/cache/pwsh.tgz -C /tmp/p
/tmp/p/pwsh -NoProfile -Command '\$e = \$null; [void][System.Management.Automation.Language.Parser]::ParseFile(\"/t/src/install.ps1\", [ref]\$null, [ref]\$e); if (\$e) { \$e; exit 1 }; \"parsed\"'"
    r=$?
  else
    r=dl; echo "PowerShell download does not match its SHA-256" > "$D/log/ps1.log"
  fi
  if [ "$r" = 0 ] && grep -q parsed "$D/log/ps1.log"; then
    echo "SKIP ps1: Windows only; install.ps1 parses (pwsh 7.7); release.yml runs it on windows-2025"
  else
    bad=1; echo "FAIL ps1 (rc=$r): $(tail -3 "$D/log/ps1.log" | tr '\n' ' ' | cut -c1-240)"
  fi
fi

if want docker; then
  { timeout 1800 "$ctr" build -q -t crimlab:tabs "$D/src" &&
    "$ctr" run --rm crimlab:tabs crimlab orrisk --or 2 --rate 0.1 --share 0.3; } > "$D/log/docker.log" 2>&1
  res docker $?
fi

if want brew; then
  s=$(sha256sum "$D/rel/download/v$V/crimlab-$V.tar.gz" | cut -d' ' -f1)
  sed -e "s|https://github.com/rootcoder007/crimlab/releases/download/v0.0.0/crimlab-0.0.0.tar.gz|$B/download/v$V/crimlab-$V.tar.gz|" \
      -e "s|sha256 \"0\{64\}\"|sha256 \"$s\"|" "$top/packaging/brew/crimlab.rb" > "$D/art/crimlab.rb"
  # the tap rootcoder007/crimlab, made locally
  run brew docker.io/homebrew/brew:latest "export HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_ANALYTICS=1 HOMEBREW_NO_INSTALL_FROM_API=1
t=\$(brew --repository)/Library/Taps/rootcoder007/homebrew-crimlab
mkdir -p \$t/Formula; cp /t/art/crimlab.rb \$t/Formula/
brew install rootcoder007/crimlab/crimlab
brew test rootcoder007/crimlab/crimlab
crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res brew $?
fi

if want scoop; then
  if python3 -c 'import json,sys; m=json.load(open(sys.argv[1])); assert m["bin"] and m["architecture"]["64bit"]["url"]' \
    "$top/packaging/scoop/crimlab.json" 2>"$D/log/scoop.log"; then
    echo "SKIP scoop: Windows only; the manifest is valid JSON with url and bin"
  else
    bad=1; echo "FAIL scoop: $(tail -1 "$D/log/scoop.log")"
  fi
fi

if want pip; then
  run build-pip docker.io/library/python:3.13 "cp -R /t/src /tmp/s; mkdir -p /t/art/py
pip wheel -q --no-deps -w /t/art/py /tmp/s/py" || echo "build-pip failed"
  run pip docker.io/library/python:3.13-slim "pip install -q --no-index --find-links /t/art/py crimlab
python -c 'import crimlab; print(crimlab.orrisk(2, 0.1, 0.3))'"
  res pip $?
fi

if want npm; then
  if get "$WASI" "$WASISUM" wasi.tgz; then
    run build-npm docker.io/library/node:24 "mkdir /tmp/w; tar -xzf /t/cache/wasi.tgz -C /tmp/w --strip-components=1
cp -R /t/src /tmp/s; cd /tmp/s/bindings/npm; mkdir -p /t/art/npm
WASI_SDK_PATH=/tmp/w sh build.sh >/dev/null; npm pack -q --pack-destination /t/art/npm" || echo "build-npm failed"
    run npm docker.io/library/node:24-slim "cd /tmp; npm install -s /t/art/npm/crimlab-*.tgz
npx crimlab orrisk --or 2 --rate 0.1 --share 0.3
echo 'import c from \"crimlab\"; console.log(JSON.stringify(c.orrisk({or: 2, rate: 0.1, share: 0.3})));' > t.mjs; node t.mjs"
    res npm $?
  else
    bad=1; echo "FAIL npm: the wasi-sdk download does not match its SHA-256"
  fi
fi

if want cargo; then
  run build-cargo docker.io/library/rust:1 "cp -R /t/src /tmp/s; cd /tmp/s/bindings/rust; sh tools/vendor.sh
cargo package -q --allow-dirty --no-verify; mkdir -p /t/art/cargo; cp target/package/*.crate /t/art/cargo/" || echo "build-cargo failed"
  run cargo docker.io/library/rust:1 "cd /tmp; mkdir c; tar -xzf /t/art/cargo/crimlab-*.crate -C c
cargo new -q t; cd t; cargo add -q crimlab --path \$(ls -d /tmp/c/crimlab-*)
echo 'fn main() { println!(\"{}\", crimlab::orrisk(r#\"{\"or\": 2, \"rate\": 0.1, \"share\": 0.3}\"#).unwrap()); }' > src/main.rs
cargo run -q"
  res cargo $?
fi

if want go; then
  # a file module proxy holding the module zip of bindings/go, as the Go
  # proxy would serve it after the tag bindings/go/v$V
  python3 -I - "$D/src/bindings/go" "$D/art/goproxy" "$V" <<'PY' > "$D/log/build-go.log" 2>&1
import os, sys, zipfile
src, out, v = sys.argv[1], sys.argv[2], "v" + sys.argv[3]
mod = "github.com/rootcoder007/crimlab/bindings/go"
d = os.path.join(out, mod, "@v")
os.makedirs(d, exist_ok=True)
with zipfile.ZipFile(os.path.join(d, v + ".zip"), "w", zipfile.ZIP_DEFLATED) as z:
    for r, _, fs in os.walk(src):
        for f in sorted(fs):
            p = os.path.join(r, f)
            z.write(p, mod + "@" + v + "/" + os.path.relpath(p, src))
open(os.path.join(d, v + ".mod"), "w").write(open(os.path.join(src, "go.mod")).read())
open(os.path.join(d, v + ".info"), "w").write('{"Version":"%s","Time":"2026-01-01T00:00:00Z"}' % v)
open(os.path.join(d, "list"), "w").write(v + "\n")
PY
  run go docker.io/library/golang:1 "export GOPROXY=file:///t/art/goproxy GOSUMDB=off GOFLAGS=-mod=mod
cd /tmp; mkdir t; cd t; go mod init t >/dev/null 2>&1
go get github.com/rootcoder007/crimlab/bindings/go@v$V
cat > main.go <<'G'
package main

import (
	\"fmt\"

	crimlab \"github.com/rootcoder007/crimlab/bindings/go\"
)

func main() {
	r, err := crimlab.Orrisk(\`{\"or\": 2, \"rate\": 0.1, \"share\": 0.3}\`)
	if err != nil {
		panic(err)
	}
	fmt.Println(r)
}
G
go run ."
  res go $?
fi

if want maven; then
  run build-maven docker.io/library/maven:3-eclipse-temurin-25 "apt-get update -qq; apt-get install -y -qq g++ make >/dev/null
cp -R /t/src /tmp/s; cd /tmp/s/bindings/java; JOBS=4 sh build.sh bundled >/dev/null
d=/t/art/m/io/github/rootcoder007/crimlab/$V; mkdir -p \$d
cp build/crimlab-bundled.jar \$d/crimlab-$V.jar; sed 's/@VERSION@/$V/' pom.xml.in > \$d/crimlab-$V.pom" || echo "build-maven failed"
  run maven docker.io/library/maven:3-eclipse-temurin-25 "cd /tmp; mkdir -p t/src/main/java; cd t
cat > pom.xml <<'M'
<project><modelVersion>4.0.0</modelVersion><groupId>t</groupId><artifactId>t</artifactId><version>1</version>
<properties><maven.compiler.release>22</maven.compiler.release></properties>
<repositories><repository><id>files</id><url>file:///t/art/m</url></repository></repositories>
<dependencies><dependency><groupId>io.github.rootcoder007</groupId><artifactId>crimlab</artifactId><version>$V</version></dependency></dependencies></project>
M
echo 'public class T { public static void main(String[] a) { System.out.println(io.github.rootcoder007.crimlab.Crimlab.jrun(\"orrisk\", \"{\\\"or\\\": 2, \\\"rate\\\": 0.1, \\\"share\\\": 0.3}\")); } }' > src/main/java/T.java
mvn -q -B dependency:build-classpath -Dmdep.outputFile=cp.txt; mvn -q -B compile
java --enable-native-access=ALL-UNNAMED -cp target/classes:\$(cat cp.txt) T"
  res maven $?
fi

if want r; then
  run build-r docker.io/library/r-base:latest "cp -R /t/src /tmp/s; cd /tmp
R CMD build --no-build-vignettes --no-manual s >/dev/null
mkdir -p /t/art/r/src/contrib; cp crimlab_*.tar.gz /t/art/r/src/contrib/
Rscript -e 'tools::write_PACKAGES(\"/t/art/r/src/contrib\", type = \"source\")'" || echo "build-r failed"
  run r docker.io/library/r-base:latest "Rscript -e 'install.packages(\"crimlab\", repos = c(\"file:///t/art/r\", \"https://cloud.r-project.org\"), quiet = TRUE); library(crimlab); print(orrisk(2, 0.1, 0.3), digits = 7)'"
  res r $?
fi

if want apt; then
  run build-deb docker.io/library/debian:trixie "apt-get update -qq; apt-get install -y -qq git build-essential >/dev/null
git config --global --add safe.directory '*'
cp -R /t/src /tmp/s; cd /tmp/s; mkdir /tmp/bd; cp -R packaging/debian /tmp/bd/; apt-get build-dep -y -qq /tmp/bd >/dev/null
sh packaging/deb.sh /t/art/deb >/dev/null" || echo "build-deb failed"
  run apt docker.io/library/debian:trixie-slim "apt-get update -qq; cd /t/art/deb; apt-get install -y -qq ./*.deb >/dev/null
crimlab -V; crimlab selftest >/dev/null; test -f /usr/include/crimlab/jrun.h
crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res apt $?
fi

if want dnf; then
  run build-rpm docker.io/library/fedora:44 "dnf -y -q install rpm-build gcc-c++ make git gawk >/dev/null
cp -R /t/src /tmp/s; cd /tmp/s; make -f .copr/Makefile srpm outdir=/tmp/o >/dev/null
rpmbuild --rebuild --define '_rpmdir /t/art/rpm' /tmp/o/*.src.rpm >/dev/null" || echo "build-rpm failed"
  run dnf docker.io/library/fedora:44 "dnf -y -q install \$(ls /t/art/rpm/*/*.rpm | grep -v debug) >/dev/null
crimlab -V; crimlab selftest >/dev/null; test -f /usr/include/crimlab/jrun.h
crimlab orrisk --or 2 --rate 0.1 --share 0.3"
  res dnf $?
fi

exit "$bad"
