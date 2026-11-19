// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the kernels
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "svla.cpp"
#include "svt.cpp"
#include "svy.cpp"
#include "swt.cpp"
#include "syn.cpp"
#include "t_aor.cpp"
#include "t_apc.cpp"
#include "t_apol.cpp"
#include "t_apx.cpp"
#include "t_buf.cpp"
#include "t_cnc.cpp"
#include "t_cnt.cpp"
#include "t_corr.cpp"
#include "t_dfit.cpp"
#include "t_dkf.cpp"
