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
