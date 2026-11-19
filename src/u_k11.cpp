// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "t_sev.cpp"
#include "t_sk.cpp"
#include "t_skch.cpp"
#include "t_sp2.cpp"
#include "t_spc.cpp"
#include "t_spi.cpp"
#include "t_svg.cpp"
#include "t_svt.cpp"
#include "t_tab.cpp"
#include "t_thr.cpp"
#include "t_tox.cpp"
#include "t_tph.cpp"
#include "t_trn.cpp"
#include "t_wdq.cpp"
#include "t_xpl.cpp"
