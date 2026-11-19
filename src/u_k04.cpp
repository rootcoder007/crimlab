// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "kmr.cpp"
#include "knx.cpp"
#include "krg.cpp"
#include "kscn.cpp"
#include "kutil.cpp"
#include "lab.cpp"
#include "lcr.cpp"
#include "learn.cpp"
#include "lnk.cpp"
#include "lnr.cpp"
#include "lsa.cpp"
#include "lws.cpp"
#include "map.cpp"
#include "mnd.cpp"
#include "mtr.cpp"
#include "mvt.cpp"
#include "net.cpp"
#include "nnd.cpp"
#include "oax.cpp"
#include "ovr.cpp"
