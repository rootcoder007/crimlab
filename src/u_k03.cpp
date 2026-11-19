// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "geou.cpp"
#include "gla.cpp"
#include "gmx.cpp"
#include "gol.cpp"
#include "gst.cpp"
#include "gsw.cpp"
#include "gtw.cpp"
#include "hot.cpp"
#include "hsa.cpp"
#include "hsp.cpp"
#include "ibn.cpp"
#include "igc.cpp"
#include "igv.cpp"
#include "jdg.cpp"
#include "jtc.cpp"
