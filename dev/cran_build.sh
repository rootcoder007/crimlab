#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Builds the CRAN tarball of crimlab: a copy of the tracked tree without the
# animated README demo (CRAN takes no such file), with the guides of
# vignettes/ combined into five (dev/cran_merge_vignettes.py), built with
# NOT_CRAN set so the five are rendered in full into the tarball; CRAN's own
# check then re-builds only their first chunk. The size is printed at the end.
#   dev/cran_build.sh [output directory, default /tmp/cran]
set -eu
src=$(cd "$(dirname "$0")/.." && pwd)
pkg=$(sed -n 's/^Package: //p' "$src/DESCRIPTION")
out=${1:-/tmp/cran}
rm -rf "${out:?}/${pkg:?}"
mkdir -p "$out/$pkg"
# only what git tracks: a stray log or object file never reaches the tarball
git -C "$src" archive HEAD | tar -xf - -C "$out/$pkg"
cd "$out/$pkg"
rm -f man/figures/demo.svg
python3 - <<'PY'
s = open("README.md", encoding="utf-8").read().split("\n")
keep = [l for l in s if "man/figures/demo.svg" not in l]
assert len(keep) == len(s) - 1, "README.md: the demo line was not found once"
open("README.md", "w", encoding="utf-8").write("\n".join(keep))
PY
python3 "$src/dev/cran_merge_vignettes.py"
cd "$out"
NOT_CRAN=true R CMD build "$pkg" > "$out/build_$pkg.log" 2>&1 || { tail -20 "$out/build_$pkg.log" >&2; exit 1; }
if tar tzf "$out"/${pkg}_*.tar.gz | grep -E '[.](o|so|dll|gcno|gcda|gcov)$|demo[.]svg$'; then
  echo "cran_build.sh: build artifacts or the demo in the tarball (above)" >&2
  exit 1
fi
ls -la "$out"/${pkg}_*.tar.gz | awk '{printf "%s %.2f MiB\n", $9, $5/1048576}'
tar tzf "$out"/${pkg}_*.tar.gz | grep -c 'inst/doc/.*[.]html$' | sed 's/^/vignettes in the tarball: /'
