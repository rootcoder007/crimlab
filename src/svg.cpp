// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab: a line chart as a standalone SVG document, the picture of the
// plots of every binding (tipplt, eceplt). No allocation: the document is
// written into the caller's buffer, or only measured. Numbers are written
// by integer arithmetic, never by printf, so the decimal point is "." in
// every locale and the bytes are the same on every platform.
#include <cmath>
#include <cstring>

#include "crimlab.h"

namespace {

// The canvas and the plot box, in pixels (documented in crimlab.h).
constexpr int W = 640, H = 400, L = 64, R = 16, T = 40, B = 52;

// Output: at most cap bytes into out (one kept for the NUL); n counts every
// byte, written or not.
struct Out {
  char *s;
  long cap;
  long n;
};

void put(Out &o, const char *t, long m) {
  for (long i = 0; i < m; ++i, ++o.n) {
    if (o.n < o.cap - 1) o.s[o.n] = t[i];
  }
}

void put(Out &o, const char *t) { put(o, t, static_cast<long>(std::strlen(t))); }

// XML text: & < > " escaped.
void text(Out &o, const char *t) {
  for (; t != nullptr && *t; ++t) {
    switch (*t) {
      case '&': put(o, "&amp;"); break;
      case '<': put(o, "&lt;"); break;
      case '>': put(o, "&gt;"); break;
      case '"': put(o, "&quot;"); break;
      default: put(o, t, 1);
    }
  }
}

// A whole number.
void whole(Out &o, long long v) {
  char b[24];
  int i = 0;
  const bool neg = v < 0;
  unsigned long long u = neg ? 0ULL - static_cast<unsigned long long>(v) : static_cast<unsigned long long>(v);
  do {
    b[i++] = static_cast<char>('0' + u % 10);
    u /= 10;
  } while (u > 0);
  if (neg) put(o, "-");
  while (i > 0) put(o, &b[--i], 1);
}

// v with d decimals (0 <= d <= 15), |v| * 10^d below 1e15; no "-0".
void fixed(Out &o, double v, int d) {
  double p = 1;
  for (int i = 0; i < d; ++i) p *= 10;
  const long long u = std::llround(v * p);
  const long long q = static_cast<long long>(p);
  if (u < 0) put(o, "-");
  const unsigned long long a = u < 0 ? 0ULL - static_cast<unsigned long long>(u) : static_cast<unsigned long long>(u);
  whole(o, static_cast<long long>(a / static_cast<unsigned long long>(q)));
  if (d > 0) {
    put(o, ".");
    char b[16];
    unsigned long long f = a % static_cast<unsigned long long>(q);
    for (int i = d - 1; i >= 0; --i) {
      b[i] = static_cast<char>('0' + f % 10);
      f /= 10;
    }
    put(o, b, d);
  }
}

// A tick label: fixed with d decimals while it fits, otherwise m.mme+E.
void label(Out &o, double v, int d) {
  if (v == 0 || (std::fabs(v) < 1e15 && std::fabs(v) * std::pow(10.0, d) < 1e15)) {
    fixed(o, v == 0 ? 0.0 : v, d);
    return;
  }
  // a tick is 1, 2 or 5 times a power of 10, so its mantissa is never
  // near 10; the 1e-9 only keeps log10's rounding from taking a power of
  // 10 for the one below it
  const int e = static_cast<int>(std::floor(std::log10(std::fabs(v)) + 1e-9));
  fixed(o, v / std::pow(10.0, e), 2);
  put(o, e < 0 ? "e-" : "e+");
  whole(o, e < 0 ? -e : e);
}

// A pixel coordinate, two decimals.
void px(Out &o, double v) { fixed(o, v, 2); }

// Ticks over [a, b]: a step of 1, 2 or 5 times a power of 10, about five
// of them; i0 and i1 the first and last multiples of the step in range, d
// the decimals a label needs.
void ticks(double a, double b, double *step, long long *i0, long long *i1, int *d) {
  const double raw = (b - a) / 5;
  const double e = std::floor(std::log10(raw));
  const double mag = std::pow(10.0, e);
  const double f = raw / mag;
  const double m = f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10;
  *step = m * mag;
  *i0 = static_cast<long long>(std::ceil(a / *step - 1e-9));
  *i1 = static_cast<long long>(std::floor(b / *step + 1e-9));
  const int k = static_cast<int>(m == 10 ? e + 1 : e);
  *d = k < 0 ? (-k > 15 ? 15 : -k) : 0;
}

// The data range [lo, hi] padded 4% each way (R's default); a single
// value is widened by 1 either way first.
void range(double *lo, double *hi) {
  if (*hi - *lo <= 0) {
    *lo -= 1;
    *hi += 1;
  }
  const double p = 0.04 * (*hi - *lo);
  *lo -= p;
  *hi += p;
}

const char *const dashes[] = {"", " stroke-dasharray=\"7 4\"", " stroke-dasharray=\"2 3\""};

}  // namespace

extern "C" int clab_svgln(int k, const int *off, const int *np, const double *x,
                          const double *y, const int *step, const int *dash,
                          const char *const *lab, const char *title, const char *xlab,
                          const char *ylab, int vline, double vx, int corner, char *out,
                          int cap, int *len) {
  if (k < 1 || k > 3 || corner < 0 || corner > 3 || cap < 0 || (cap > 0 && out == nullptr)) {
    return CLAB_EARG;
  }
  double x0 = INFINITY, x1 = -INFINITY, y0 = 0, y1 = -INFINITY;
  for (int j = 0; j < k; ++j) {
    if (np[j] < 1 || off[j] < 0 || dash[j] < 0 || dash[j] > 2) return CLAB_EARG;
    for (int i = off[j]; i < off[j] + np[j]; ++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
      if (x[i] < x0) x0 = x[i];
      if (x[i] > x1) x1 = x[i];
      if (y[i] < y0) y0 = y[i];
      if (y[i] > y1) y1 = y[i];
    }
  }
  if (vline && !std::isfinite(vx)) return CLAB_EDOM;
  if (vline && vx < x0) x0 = vx;
  if (vline && vx > x1) x1 = vx;
  range(&x0, &x1);
  range(&y0, &y1);
  const double pw = W - L - R, ph = H - T - B;
  auto mx = [&](double v) { return L + (v - x0) / (x1 - x0) * pw; };
  auto my = [&](double v) { return (H - B) - (v - y0) / (y1 - y0) * ph; };
  Out o{out, cap, 0};
  put(o, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"400\" "
         "viewBox=\"0 0 640 400\" role=\"img\">\n<title>");
  text(o, title);
  put(o, "</title>\n<style>\n"
         "svg{font-family:system-ui,-apple-system,\"Segoe UI\",Roboto,\"Helvetica Neue\","
         "Arial,sans-serif;font-size:12px}\n"
         "text{fill:#3a3a3a}\n.ax{stroke:#3a3a3a;fill:none}\n"
         ".gr{stroke:#8a8a8a;stroke-opacity:.35;fill:none}\n"
         ".c0{stroke:#2f72b8;fill:none;stroke-width:2}\n"
         ".c1{stroke:#c9661f;fill:none;stroke-width:2}\n"
         ".c2{stroke:#4c9a55;fill:none;stroke-width:2}\n"
         "@media (prefers-color-scheme:dark){text{fill:#d6d6d6}.ax{stroke:#d6d6d6}"
         ".c0{stroke:#6fa8e6}.c1{stroke:#ec9a5c}.c2{stroke:#7cc485}}\n"
         "</style>\n<g class=\"plot\">\n");
  // grid, ticks and their labels
  for (int a = 0; a < 2; ++a) {
    double st;
    long long i0, i1;
    int d;
    if (a == 0) ticks(x0, x1, &st, &i0, &i1, &d); else ticks(y0, y1, &st, &i0, &i1, &d);
    for (long long i = i0; i <= i1; ++i) {
      const double v = static_cast<double>(i) * st;
      put(o, "<path class=\"gr\" d=\"M");
      if (a == 0) {
        px(o, mx(v)); put(o, " "); px(o, T); put(o, "V"); px(o, H - B);
      } else {
        px(o, L); put(o, " "); px(o, my(v)); put(o, "H"); px(o, W - R);
      }
      put(o, "\"/>\n<text ");
      if (a == 0) {
        put(o, "x=\""); px(o, mx(v)); put(o, "\" y=\""); px(o, H - B + 16);
        put(o, "\" text-anchor=\"middle\">");
      } else {
        put(o, "x=\""); px(o, L - 6); put(o, "\" y=\""); px(o, my(v) + 4);
        put(o, "\" text-anchor=\"end\">");
      }
      label(o, v, d);
      put(o, "</text>\n");
    }
  }
  put(o, "<path class=\"ax\" d=\"M");
  px(o, L); put(o, " "); px(o, T); put(o, "V"); px(o, H - B); put(o, "H"); px(o, W - R);
  put(o, "V"); px(o, T); put(o, "Z\"/>\n");
  if (vline) {
    put(o, "<path class=\"ax vline\" stroke-dasharray=\"2 3\" d=\"M");
    px(o, mx(vx)); put(o, " "); px(o, T); put(o, "V"); px(o, H - B); put(o, "\"/>\n");
  }
  // the curves, one path each
  for (int j = 0; j < k; ++j) {
    put(o, "<path class=\"curve c");
    whole(o, j);
    put(o, "\"");
    put(o, dashes[dash[j]]);
    put(o, " d=\"M");
    const int a = off[j], e = off[j] + np[j];
    px(o, mx(x[a])); put(o, " "); px(o, my(y[a]));
    for (int i = a + 1; i < e; ++i) {
      if (step[j]) {
        put(o, "H"); px(o, mx(x[i])); put(o, "V"); px(o, my(y[i]));
      } else {
        put(o, "L"); px(o, mx(x[i])); put(o, " "); px(o, my(y[i]));
      }
    }
    put(o, "\"/>\n");
  }
  // legend, in the corner asked for, inside the box
  put(o, "<g class=\"legend\">\n");
  const double lx = (corner & 1) ? W - R - 170 : L + 12;
  const double ly0 = (corner & 2) ? H - B - 16 - 18 * (k - 1) : T + 16;
  for (int j = 0; j < k; ++j) {
    const double ly = ly0 + 18 * j;
    put(o, "<line class=\"c");
    whole(o, j);
    put(o, "\"");
    put(o, dashes[dash[j]]);
    put(o, " x1=\""); px(o, lx); put(o, "\" y1=\""); px(o, ly);
    put(o, "\" x2=\""); px(o, lx + 30); put(o, "\" y2=\""); px(o, ly);
    put(o, "\"/>\n<text x=\""); px(o, lx + 36); put(o, "\" y=\""); px(o, ly + 4);
    put(o, "\">");
    text(o, lab[j]);
    put(o, "</text>\n");
  }
  put(o, "</g>\n</g>\n<text class=\"title\" x=\"");
  px(o, L); put(o, "\" y=\"24\" font-size=\"15\" font-weight=\"600\">");
  text(o, title);
  put(o, "</text>\n<text class=\"xlab\" x=\"");
  px(o, L + pw / 2); put(o, "\" y=\""); px(o, H - 12);
  put(o, "\" text-anchor=\"middle\">");
  text(o, xlab);
  put(o, "</text>\n<text class=\"ylab\" transform=\"translate(16 ");
  px(o, T + ph / 2); put(o, ") rotate(-90)\" text-anchor=\"middle\">");
  text(o, ylab);
  put(o, "</text>\n</svg>\n");
  if (o.n > 2147483646L) return CLAB_EARG;
  *len = static_cast<int>(o.n);
  if (cap > 0) {
    if (o.n >= cap) return CLAB_EARG;
    out[o.n] = '\0';
  }
  return CLAB_OK;
}
