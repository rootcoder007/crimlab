// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "t_hwk.cpp"
#include "t_hyq.cpp"
#include "t_ibn.cpp"
#include "t_igh.cpp"
#include "t_jdg.cpp"
#include "t_jtc.cpp"
#include "t_kmr.cpp"
#include "t_krg.cpp"
#include "t_lab.cpp"
#include "t_lcr.cpp"
#include "t_lnk.cpp"
#include "t_lnr.cpp"
#include "t_lws.cpp"
#include "t_map.cpp"
#include "t_mhk.cpp"
#include "t_mnd.cpp"
#include "t_mtr.cpp"
#include "t_net.cpp"
