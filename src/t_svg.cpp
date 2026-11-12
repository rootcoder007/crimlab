// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The SVG line chart (svg.cpp, clab_svgln): two passes agree, the bytes
// repeat, text is escaped, the data map to the plot box, and bad input is
// refused. Built only with the self-test (selftest.cpp).
#include <cstring>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

char buf[16384];

// The number after key (" d=\"M" and the like) in s, read by hand.
double after(const char *s, const char *key) {
  const char *p = std::strstr(s, key);
  if (p == nullptr) return -1;
  p += std::strlen(key);
  double v = 0, f = 0.1;
  bool frac = false;
  for (; (*p >= '0' && *p <= '9') || *p == '.'; ++p) {
    if (*p == '.') frac = true;
    else if (!frac) v = v * 10 + (*p - '0');
    else { v += (*p - '0') * f; f /= 10; }
  }
  return v;
}

}  // namespace

void t_svg(tk &t, int) {
  const double x[] = {-1, 0, 2, -0.5, 1.5};
  const double y[] = {1, 0.5, 0.25, 0.75, 0};
  const int off[] = {0, 3}, np[] = {3, 2}, step[] = {1, 0}, dash[] = {0, 1};
  const char *const lab[] = {"a & b", "<c>"};
  const char *ti = "t \"q\"";
  int n = -1, m = -1;
  tk_ieq(t, "svg1",
         clab_svgln(2, off, np, x, y, step, dash, lab, ti, "x", "y", 1, 0, 1, nullptr, 0, &n),
         CLAB_OK);
  tk_true(t, "svg2", n > 1000 && n < static_cast<int>(sizeof buf));
  tk_ieq(t, "svg3",
         clab_svgln(2, off, np, x, y, step, dash, lab, ti, "x", "y", 1, 0, 1, buf, n + 1, &m),
         CLAB_OK);
  tk_ieq(t, "svg4", m, n);
  tk_ieq(t, "svg5", static_cast<int>(std::strlen(buf)), n);
  tk_true(t, "svg6", std::strstr(buf, "a &amp; b") && std::strstr(buf, "&lt;c&gt;") &&
                         std::strstr(buf, "t &quot;q&quot;"));
  tk_true(t, "svg7", std::strncmp(buf + n - 7, "</svg>\n", 7) == 0);
  // x spans -1 .. 2 padded 4%: -1.12 .. 2.12; y spans 0 .. 1: -0.04 .. 1.04.
  // The first curve starts at (-1, 1): x 64 + 0.12 / 3.24 * 560 = 84.74.
  const char *c0 = std::strstr(buf, "class=\"curve c0\"");
  tk_true(t, "svg8", c0 != nullptr);
  tk_near(t, "svg9", after(c0, " d=\"M"), 64 + 0.12 / 3.24 * 560, 0.006);
  tk_near(t, "svg10", after(c0, " d=\"M"), 84.74, 0.006);
  // the bytes repeat
  static char again[16384];
  clab_svgln(2, off, np, x, y, step, dash, lab, ti, "x", "y", 1, 0, 1, again, n + 1, &m);
  tk_true(t, "svg11", std::memcmp(buf, again, static_cast<std::size_t>(n) + 1) == 0);
  // one point: widened by 1 either way
  const double one = 3;
  const int o1[] = {0}, n1[] = {1};
  tk_ieq(t, "svg12",
         clab_svgln(1, o1, n1, &one, &one, step, dash, lab, "", "", "", 0, 0, 0, nullptr, 0, &m),
         CLAB_OK);
  // the legend in the corner asked for: bottom left starts at x 76
  clab_svgln(2, off, np, x, y, step, dash, lab, ti, "x", "y", 1, 0, 2, buf, n + 1, &m);
  const char *lg = std::strstr(buf, "class=\"legend\"");
  tk_near(t, "svg13", after(lg, " x1=\""), 76, 0.006);
  // ticks too large for fixed notation are labelled m.mme+E
  const double big[] = {0, 4e20}, by[] = {0, 1};
  const int n2b[] = {2};
  clab_svgln(1, o1, n2b, big, by, step, dash, lab, "", "", "", 0, 0, 0, buf,
             static_cast<int>(sizeof buf), &m);
  tk_true(t, "svg20", std::strstr(buf, "e+20<") != nullptr);
  // refused: a cap too small, no curve, an empty one, a bad corner,
  // non-finite values
  tk_ieq(t, "svg14",
         clab_svgln(2, off, np, x, y, step, dash, lab, "", "", "", 0, 0, 0, buf, 10, &m),
         CLAB_EARG);
  tk_ieq(t, "svg15",
         clab_svgln(0, off, np, x, y, step, dash, lab, "", "", "", 0, 0, 0, nullptr, 0, &m),
         CLAB_EARG);
  const int n0[] = {0};
  tk_ieq(t, "svg16",
         clab_svgln(1, o1, n0, x, y, step, dash, lab, "", "", "", 0, 0, 0, nullptr, 0, &m),
         CLAB_EARG);
  tk_ieq(t, "svg17",
         clab_svgln(1, o1, n1, x, y, step, dash, lab, "", "", "", 0, 0, 4, nullptr, 0, &m),
         CLAB_EARG);
  const double inf[] = {1, std::numeric_limits<double>::infinity()};
  const int n2[] = {2};
  tk_ieq(t, "svg18",
         clab_svgln(1, o1, n2, x, inf, step, dash, lab, "", "", "", 0, 0, 0, nullptr, 0, &m),
         CLAB_EDOM);
  tk_ieq(t, "svg19",
         clab_svgln(1, o1, n2, x, y, step, dash, lab, "", "", "", 1,
                    std::numeric_limits<double>::quiet_NaN(), 0, nullptr, 0, &m),
         CLAB_EDOM);
}
