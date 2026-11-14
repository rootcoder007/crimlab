#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Prints the preset table of the R package (inst/presets/presets.tsv), the
# rows that cli/pre.inc holds, one {9 fields} per line after the header:
#   sh cli/pre.sh > cli/pre.inc
# CI checks that pre.inc is up to date (pkg.yml).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
LC_ALL=C awk -F "\t" "NR > 1 {
  sub(/\r\$/, \"\")
  o = \"    {\"
  for (i = 1; i <= 9; i++) {
    v = i <= NF ? \$i : \"\"
    gsub(/\\\\/, \"\\\\\\\\\", v)
    gsub(/\"/, \"\\\\\\\"\", v)
    o = o (i > 1 ? \", \" : \"\") \"\\\"\" v \"\\\"\"
  }
  print o \"},\"
}" "${1:-$here/../inst/presets/presets.tsv}"
