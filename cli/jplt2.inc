// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// More plots of R/plotx.R for clab_jrun(), as in jplt.inc: nothing is
// drawn, and each returns the rows the R function returns invisibly, with
// the checks of R in their order and words. The densities of plotrdg()
// and plotvln() are clab_dens() (src/dens.cpp), stats::density() with R's
// bandwidth rules, binning, FFT and interpolation: equal to R's to
// rounding, not bit for bit. The attribute `summary` of plotvln() is not
// part of the result.

namespace plt {

// .pbw(): the rule (0 for a number) and the bandwidth.
inline bool pbw(Ctx &c, int &rule, double &bw) {
  rule = 1;
  bw = 0.0;
  if (!c.has("bw")) return true;
  const Val &b = c["bw"];
  static const char *const rs[] = {"nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "sj-dpi"};
  static const int rn[] = {1, 2, 3, 4, 5, 5, 6};
  if (b.t == Val::STR) {
    for (std::size_t i = 0; i < 7; ++i) {
      if (b.s == rs[i]) {
        rule = rn[i];
        return true;
      }
    }
  } else if (ax::numv(b)) {
    const std::vector<double> v = vals(b);
    if (v.size() == 1 && std::isfinite(v[0]) && v[0] > 0) {
      rule = 0;
      bw = v[0];
      return true;
    }
  }
  return c.fail("`bw` must be one positive number or one of nrd0, nrd, ucv, bcv, sj, sj-ste, sj-dpi.");
}

// stats::density(v, bw = bw): its 512 points and densities, or R's
// message (without the name of this function, as R gives it).
inline bool dens(Ctx &c, const std::vector<double> &v, int rule, double bw, std::vector<double> &x,
                 std::vector<double> &y) {
  x.assign(512, 0.0);
  y.assign(512, 0.0);
  double st[2] = {0.0, 0.0};
  int why = 0;
  const int rc = clab_dens(v.data(), nullptr, static_cast<std::int64_t>(v.size()), rule, bw, 1.0, 0, 512, 3.0, 0, 0.0,
                           0, 0.0, x.data(), y.data(), st, &why);
  if (rc == CLAB_EDOM) {
    static const char *const m[] = {"",
                                    "need at least 2 points to select a bandwidth automatically",
                                    "data are constant in bandwidth calculation",
                                    "sample is too sparse to find TD",
                                    "sample is too sparse to find alph2",
                                    "no solution in the specified range of bandwidths",
                                    "non-finite 'bw'",
                                    "'bw' is not positive.",
                                    "non-finite 'from'"};
    c.err = m[why >= 0 && why <= 8 ? why : 0];
    return false;
  }
  return kchk(c, rc);
}

// The values of x in group l.
inline std::vector<double> ingrp(const Xg &d, const std::string &l) {
  std::vector<double> v;
  for (std::size_t i = 0; i < d.x.size(); ++i) {
    if (d.g[i] == l) v.push_back(d.x[i]);
  }
  return v;
}

// The densities of the groups with two values or more, as rows.
inline bool grpdens(Ctx &c, const Xg &d, int rule, double bw, Frame &f, bool &any) {
  ax::Cols g;
  std::vector<double> xs, ys;
  any = false;
  for (const std::string &l : d.lv) {
    const std::vector<double> v = ingrp(d, l);
    if (v.size() < 2) continue;
    std::vector<double> x, y;
    CK(dens(c, v, rule, bw, x, y));
    any = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
      g.push_back(l);
      xs.push_back(x[i]);
      ys.push_back(y[i]);
    }
  }
  f.str("grp", g);
  f.num("x", xs);
  f.num("y", ys);
  return true;
}

}  // namespace plt

// plotrdg()
bool plotrdg(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  int rule = 1;
  double bw = 0.0;
  CK(plt::pbw(c, rule, bw));
  Frame f;
  bool any = false;
  CK(plt::grpdens(c, d, rule, bw, f, any));
  if (!any) return c.fail("no group has two values.");
  out = f.o;
  return true;
}

// plotvln()
bool plotvln(Ctx &c, Val &out) {
  plt::Xg d;
  CK(plt::pxg(c, d));
  int rule = 1;
  double bw = 0.0;
  CK(plt::pbw(c, rule, bw));
  Frame f;
  bool any = false;
  CK(plt::grpdens(c, d, rule, bw, f, any));
  if (!any) {
    // R sets the attribute summary on NULL
    c.err = "attempt to set an attribute on NULL";
    return false;
  }
  out = f.o;
  return true;
}

namespace plt {

// .prange() of finite values: their range, widened by 0.5 each side when
// it is one value.
inline void prange(const std::vector<double> &v, double &lo, double &hi) {
  lo = *std::min_element(v.begin(), v.end());
  hi = *std::max_element(v.begin(), v.end());
  if (hi - lo == 0) {
    lo -= 0.5;
    hi += 0.5;
  }
}
// R's x %% 2 of a whole number: 0 or 1.
inline double mod2(double r) { return r - 2 * std::floor(r / 2); }

}  // namespace plt

// plothex(): the hexagons in order of row, then column, as R orders them.
bool plothex(Ctx &c, Val &out) {
  const Val &xo = c["x"], &yo = c["y"];
  if (!ax::numv(xo)) return c.fail("`x` must be a numeric vector.");
  const std::vector<double> x0 = plt::vals(xo);
  std::size_t nf = 0;
  for (double v : x0) nf += std::isfinite(v) ? 1 : 0;
  if (nf < 2) return c.fail("`x` needs at least two finite values.");
  const std::vector<double> y0 = ax::numv(yo) ? plt::vals(yo) : std::vector<double>();
  if (!ax::numv(yo) || y0.size() != x0.size()) return c.fail("`y` must be numbers, one per value of `x`.");
  std::vector<double> x, y;
  for (std::size_t i = 0; i < x0.size(); ++i) {
    if (!std::isfinite(x0[i]) || !std::isfinite(y0[i])) continue;
    x.push_back(x0[i]);
    y.push_back(y0[i]);
  }
  if (x.size() < 2) return c.fail("needs at least two points with both values.");
  double bins = 30;
  if (c.has("bins")) {
    const std::vector<double> b = ax::numv(c["bins"]) ? plt::vals(c["bins"]) : std::vector<double>();
    const bool ok = b.size() == 1 && std::isfinite(b[0]) && b[0] >= 2 && b[0] <= 500 && b[0] == std::round(b[0]);
    if (!ok) return c.fail("`bins` must be one whole number from 2 to 500.");
    bins = b[0];
  }
  double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  plt::prange(x, x1, x2);
  plt::prange(y, y1, y2);
  const double h = std::sqrt(3.0) / 2;
  std::map<std::pair<double, double>, int> cnt;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double sx = (x[i] - x1) / (x2 - x1) * bins, sy = (y[i] - y1) / (y2 - y1) * bins;
    const double r1 = std::nearbyint(sy / h), c1 = std::nearbyint(sx - plt::mod2(r1) / 2);
    double r2 = std::floor(sy / h);
    if (r2 == r1) r2 += 1;
    const double c2 = std::nearbyint(sx - plt::mod2(r2) / 2);
    const double d1 = (sx - c1 - plt::mod2(r1) / 2) * (sx - c1 - plt::mod2(r1) / 2) + (sy - r1 * h) * (sy - r1 * h);
    const double d2 = (sx - c2 - plt::mod2(r2) / 2) * (sx - c2 - plt::mod2(r2) / 2) + (sy - r2 * h) * (sy - r2 * h);
    ++cnt[d1 <= d2 ? std::make_pair(r1, c1) : std::make_pair(r2, c2)];
  }
  std::vector<double> hx, hy;
  std::vector<int> n;
  for (const auto &e : cnt) {
    const double rr = e.first.first, cc = e.first.second;
    hx.push_back(x1 + (cc + plt::mod2(rr) / 2) / bins * (x2 - x1));
    hy.push_back(y1 + rr * h / bins * (y2 - y1));
    n.push_back(e.second);
  }
  Frame f;
  f.num("x", hx);
  f.num("y", hy);
  f.ints("count", n.data(), n.size());
  out = f.o;
  return true;
}
