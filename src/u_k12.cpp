// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "t_xw.cpp"
#include "tab.cpp"
#include "tfr.cpp"
#include "thr.cpp"
#include "tlv.cpp"
#include "tox.cpp"
#include "tph.cpp"
#include "tps.cpp"
#include "trj.cpp"
#include "trn.cpp"
#include "trs.cpp"
#include "tsd.cpp"
#include "txt.cpp"
#include "vrfp.cpp"
#include "wcb.cpp"
#include "wdq.cpp"
#include "xpl.cpp"
#include "xw.cpp"
#include "apx.cpp"
