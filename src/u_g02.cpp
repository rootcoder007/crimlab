// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One translation unit of the R package's build (src/Makevars): the R glue
// below, compiled together so that the headers are read once for the
// group and not once per file. Made by dev/unity.py; a new source file
// is added to a group there. The header below marks the file as part of
// the R build for the other builds, which take each source by itself.
#include <R_ext/Visibility.h>

#include "bhsa.cpp"
#include "bibn.cpp"
#include "bigh.cpp"
#include "bind.cpp"
#include "bjdg.cpp"
#include "blcr.cpp"
#include "blnk.cpp"
#include "blnr.cpp"
#include "bmap.cpp"
#include "bmat.cpp"
#include "bmnd.cpp"
#include "bmtr.cpp"
#include "bod.cpp"
#include "bplp.cpp"
#include "bplr.cpp"
#include "bprh.cpp"
#include "bprim.cpp"
#include "bpsc.cpp"
#include "bpub.cpp"
#include "bpwr.cpp"
