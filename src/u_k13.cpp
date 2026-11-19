// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "dscd.cpp"
#include "dwn.cpp"
#include "geo2.cpp"
#include "igf.cpp"
#include "igh.cpp"
#include "t_bnf.cpp"
#include "t_cgr.cpp"
