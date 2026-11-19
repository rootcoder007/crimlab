// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "ped.cpp"
#include "plf.cpp"
#include "plp.cpp"
#include "plr.cpp"
#include "ppx.cpp"
#include "prh.cpp"
#include "prim.cpp"
#include "psc.cpp"
#include "pub.cpp"
#include "pwa.cpp"
#include "pwr.cpp"
#include "pwrd.cpp"
#include "pwrmac.cpp"
#include "pwrwin.cpp"
#include "pws.cpp"
#include "qn.cpp"
#include "rct.cpp"
#include "rdd.cpp"
#include "rf.cpp"
#include "rnd.cpp"
#include "rpk.cpp"
