#!/bin/sh
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Prints the proof table of the R binding (src/bind.cpp), the rows that
# cli/thm.inc holds:   sh cli/thm.sh > cli/thm.inc
# CI checks that thm.inc is up to date (cli.yml).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
sed -n '/^const Thm thms\[\] = {/,/^};/p' "${1:-$here/../src/bind.cpp}" | grep '^ *{"'
