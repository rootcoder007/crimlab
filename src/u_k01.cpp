// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "ahl.cpp"
#include "aor.cpp"
#include "apc.cpp"
#include "apol.cpp"
#include "arl.cpp"
#include "bnf.cpp"
#include "buf.cpp"
#include "cal.cpp"
#include "cba.cpp"
#include "cgeo.cpp"
#include "cgr.cpp"
#include "chm.cpp"
#include "cmps.cpp"
#include "cnc.cpp"
#include "cnt.cpp"
#include "core.cpp"
#include "dag.cpp"
#include "dens.cpp"
