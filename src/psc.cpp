// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Credit for pre-sentence custody (R/psc.R): Criminal Code s. 719(3) and
// (3.1). Each day held before sentence is credited at a ratio of at most
// one day (s. 719(3)) or, where the circumstances justify it, one and a
// half days (s. 719(3.1)); the credit is subtracted from the sentence,
// which cannot fall below nothing.
#include <cmath>
#include <cstdint>

#include "crimlab.h"

extern "C" int clab_pscred(const double *days, int64_t n, const double *ratio, int64_t nr, const double *sent,
                           int64_t ns, double *credit, double *remain) {
  if (n < 1 || (nr != 1 && nr != n) || (sent != nullptr && ns != 1 && ns != n)) return CLAB_EARG;
  for (int64_t i = 0; i < n; ++i) {
    const double d = days[i], r = ratio[nr == 1 ? 0 : i];
    if (!std::isfinite(d) || d < 0 || !std::isfinite(r) || r < 0 || r > 1.5) return CLAB_EDOM;
    credit[i] = d * r;
    if (sent == nullptr) continue;
    const double s = sent[ns == 1 ? 0 : i];
    if (!std::isfinite(s) || s < 0) return CLAB_EDOM;
    remain[i] = s > credit[i] ? s - credit[i] : 0.0;
  }
  return CLAB_OK;
}
