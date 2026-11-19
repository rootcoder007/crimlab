// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the R glue
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "bahl.cpp"
#include "bapc.cpp"
#include "bcba.cpp"
#include "bcnc.cpp"
#include "bdfit.cpp"
#include "bdid.cpp"
#include "bdkf.cpp"
#include "bdpn.cpp"
#include "bdsc.cpp"
#include "bdscd.cpp"
#include "bdtr.cpp"
#include "bdwn.cpp"
#include "benv.cpp"
#include "bflw.cpp"
#include "bfmt.cpp"
#include "bfp.cpp"
#include "bfpio.cpp"
#include "bgol.cpp"
