// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "rpt.cpp"
#include "rpv.cpp"
#include "rsc.cpp"
#include "rtm.cpp"
#include "sar.cpp"
#include "sdb.cpp"
#include "selftest.cpp"
#include "sens.cpp"
#include "sev.cpp"
#include "sfx.cpp"
#include "sk.cpp"
#include "skch.cpp"
#include "spc.cpp"
#include "spr.cpp"
#include "spt.cpp"
#include "sreg.cpp"
#include "stg.cpp"
#include "stm.cpp"
#include "svb.cpp"
#include "svg.cpp"
