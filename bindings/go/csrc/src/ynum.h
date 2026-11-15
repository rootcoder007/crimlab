// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The numbers of a change table as the writers (R/yoyw.R) write them:
// every digit (the shortest text that reads back as the same double), a
// fixed number of decimals with commas between thousands (formatC(x,
// "f", digits, big.mark = ",")), or a signed percentage. Header only.
#ifndef CLAB_YNUM_H
#define CLAB_YNUM_H

#include <cmath>
#include <cstdio>
#include <string>

#include "dfmt.h"

namespace ynum {

inline std::string all(double x) {
  char b[40];
  const int n = dfmt::put(x, b);
  return std::string(b, static_cast<std::size_t>(n));
}

inline std::string fix(double x, int d) {
  char b[400];
  std::snprintf(b, sizeof b, "%.*f", d, x);
  std::string s = b;
  if (!std::isfinite(x)) return s;
  const std::size_t a = s[0] == '-' ? 1 : 0, e = s.find('.') == std::string::npos ? s.size() : s.find('.');
  for (std::size_t k = e; k > a + 3; k -= 3) s.insert(k - 3, ",");
  return s;
}

inline std::string pct(double x, int d, const char *suf) {
  if (std::isinf(x)) return "+Inf";
  char b[400];
  std::snprintf(b, sizeof b, "%+.*f%s", d, x, suf);
  return b;
}

}  // namespace ynum

#endif
