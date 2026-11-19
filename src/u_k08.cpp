// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "t_dml.cpp"
#include "t_dpn.cpp"
#include "t_drf.cpp"
#include "t_dsc.cpp"
#include "t_dst.cpp"
#include "t_dtr.cpp"
#include "t_dwn.cpp"
#include "t_env.cpp"
#include "t_flw.cpp"
#include "t_fp.cpp"
#include "t_fsim.cpp"
#include "t_geou.cpp"
#include "t_gmx.cpp"
#include "t_gol.cpp"
#include "t_hsa.cpp"
