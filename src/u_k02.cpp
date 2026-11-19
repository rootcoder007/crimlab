// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "did.cpp"
#include "dkf.cpp"
#include "dml.cpp"
#include "dpn.cpp"
#include "dpnd.cpp"
#include "drf.cpp"
#include "dsc.cpp"
#include "dtr.cpp"
#include "env.cpp"
#include "eof.cpp"
#include "fair.cpp"
#include "fls.cpp"
#include "flw.cpp"
#include "fp.cpp"
#include "fpdat.cpp"
#include "fsim.cpp"
#include "geo.cpp"
