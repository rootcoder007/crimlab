#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Prints the publishers' codes of the R package (inst/presets/codes.tsv),
# the rows that cli/codes.inc holds, one {4 fields} per line after the
# header:
#   sh cli/codes.sh > cli/codes.inc
# CI checks that codes.inc is up to date (pkg.yml).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
LC_ALL=C awk -F "\t" "NR > 1 {
  sub(/\r\$/, \"\")
  o = \"    {\"
  for (i = 1; i <= 4; i++) {
    v = i <= NF ? \$i : \"\"
    gsub(/\\\\/, \"\\\\\\\\\", v)
    gsub(/\"/, \"\\\\\\\"\", v)
    o = o (i > 1 ? \", \" : \"\") \"\\\"\" v \"\\\"\"
  }
  print o \"},\"
}" "${1:-$here/../inst/presets/codes.tsv}"
