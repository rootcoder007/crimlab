// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Runner for the standalone kernel tests. Built only with -DCLAB_MAIN:
//   g++ -std=c++17 -DCLAB_MAIN src/<family>.cpp src/t_<family>.cpp
//       src/t_hmain.cpp -o th && ./th        (./th bench: the timings)
// or every family at once: GATE=kern sh dev/gate.sh
#ifdef CLAB_MAIN

#include <cstdio>
#include <cstring>

#include "tk.h"

int main(int argc, char **argv) {
  const bool bench = argc > 1 && std::strcmp(argv[1], "bench") == 0;
  for (int i = 0; i < tk::ntab(); ++i) {
    const tk::ent &e = tk::tab()[i];
    if (bench) {
      if (e.bench != nullptr) e.bench();
    } else {
      e.test();
    }
  }
  if (!bench) {
    std::printf("tk: %d families, %d checks, %d failed\n", tk::ntab(),
                tk::nrun(), tk::nfail());
  }
  return tk::nfail() == 0 ? 0 : 1;
}

#endif
