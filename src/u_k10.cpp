// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "t_nnd.cpp"
#include "t_oax.cpp"
#include "t_ovr.cpp"
#include "t_plr.cpp"
#include "t_pol.cpp"
#include "t_prh.cpp"
#include "t_prim.cpp"
#include "t_pub.cpp"
#include "t_pwr.cpp"
#include "t_rat.cpp"
#include "t_rf.cpp"
#include "t_rpv.cpp"
#include "t_rtm.cpp"
#include "t_sar.cpp"
#include "t_sdb.cpp"
#include "t_sep.cpp"
