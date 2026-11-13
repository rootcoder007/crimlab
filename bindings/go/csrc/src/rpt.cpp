// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, policing: agencies that report in every year of a panel.
#include <climits>

#include "crimlab.h"

extern "C" {

int clab_rptc(const int *unit, const int *year, const int *rep, int n, int a,
              int y0, int y1, int minrep, int *flag, int *nc, int *work) {
  if (n < 0 || a < 1 || y1 < y0) return CLAB_EARG;
  if (y1 - y0 >= INT_MAX / a) return CLAB_EARG;  // a * years must fit
  const int ny = y1 - y0 + 1;
  for (int i = 0; i < n; ++i) {
    if (unit[i] < 0 || unit[i] >= a) return CLAB_EARG;
  }
  for (long long k = 0; k < static_cast<long long>(a) * ny; ++k) work[k] = 0;
  // mark (unit, year) when some record of it reaches minrep months
  for (int i = 0; i < n; ++i) {
    if (year[i] < y0 || year[i] > y1 || rep[i] < minrep) continue;
    work[static_cast<long long>(unit[i]) * ny + (year[i] - y0)] = 1;
  }
  int c = 0;
  for (int u = 0; u < a; ++u) {
    int all = 1;
    for (int y = 0; y < ny && all; ++y) all = work[static_cast<long long>(u) * ny + y];
    flag[u] = all;
    c += all;
  }
  *nc = c;
  return CLAB_OK;
}

}  // extern "C"
